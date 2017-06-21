{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, DeriveFunctor #-}
{-# LANGUAGE FlexibleContexts, MultiWayIf, TypeFamilies #-}
module ECC.Code.LDPC.GPU.CachedMult where

-- Uses the StableDiv data structure to cache multiplications.
-- Uses Accelerate to run on the GPU

import Prelude hiding ((==), (/=), (>=), (<), (>), all, map, (||), (&&), not, Num, snd)
import qualified Prelude as P

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import qualified Data.Matrix as M
import Data.Bit
import Data.Bits
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.QuasiCyclic as Q
import Debug.Trace

import qualified ECC.Code.LDPC.Fast.Encoder as E

import Data.Semigroup
import Data.Foldable (foldl')

import Data.Array.Accelerate
import Data.Array.Accelerate.IO
import Data.Array.Accelerate.LLVM.PTX

code :: Code
code = mkLDPC_Code "gpu-arraylet-cm" E.encoder decoder

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a rate maxIterations orig_lam =
  Just $ fromAcc (ldpc mLet maxIterations (toAcc sh orig_lam))
  where mLet = initMatrixlet 0 a
        sh   = Z :. U.length orig_lam

---------------------------------------------------------------------


type StableDiv
  = Exp ((,)
           Double  -- | The "worst" value for division: the closest to zero
           Double) -- | The result of the multiplications,
                   --   excluding the "worst" value

absMinMax :: Exp Double -> Exp Double -> Exp (Double, Double)
absMinMax x y =
  cond (abs x < abs y)
       (lift (x, y))
       (lift (y, x))

lit :: Exp Double -> StableDiv
lit x =
  cond (x >= 1)
       (lift (1 :: Double, x))
       (lift (x, 1 :: Double))

smult :: StableDiv -> StableDiv -> StableDiv
smult p q =
    lift (minOfMins, (b * maxOfMins * d))
    where
      minOfMins, maxOfMins :: Exp Double
      (minOfMins, maxOfMins) = unlift (absMinMax a c)

      (a, b) = unlift p
      (c, d) = unlift q

sdiv :: StableDiv -> Exp Double -> Exp Double
sdiv p c =
  cond (a == c)
       b
       (a * (b/c))
  where
    (a, b) = unlift p

type V a = Array DIM1 a
type M a = Array DIM2 a

type Arraylet a =
  (,)
    (Scalar Int)
    (V a)

arrayArraylet :: Acc (Scalar Int) -> Exp Int -> Exp Int -> (Exp DIM2 -> Exp b) -> Exp b
arrayArraylet sz0 off arrayletIx k =
  k (index2 arrayletIx ((arrayletIx + off) `mod` sz))
  where
    sz = the sz0

type Matrixlet a =
  (,,,)
    (Scalar Int)
    (M Bool)  -- | (Block) Indices with non-zero blocks
    (M Int)   -- | Offsets of arraylets
    (Array DIM3 a)

-- TODO: See if this can be made more efficient (currently traverses matrix
-- multiple times).
initMatrixlet :: Double -> Q.QuasiCyclic Integer -> Acc (Matrixlet Double)
initMatrixlet z (Q.QuasiCyclic sz m) =
  lift (unit (lift sz), indices, offsets, mat)
  where
    accM :: Acc (M Int)
    accM = lift $ fromVectors (Z :. M.nrows m :. M.ncols m) (U.convert (fmap g (M.getMatrixAsVector m)))
      where
        -- The number must be a power of two, because there is only one bit set.
        g :: Integer -> Int
        g x | x `testBit` 0 = 0
            | x P.== 0      = 0 -- Should only happen when `g` is originally called with 0
            | otherwise = 1 + g (x `shiftR` 1)

    indices :: Acc (M Bool)
    indices = map (/= 0) accM

    offsets :: Acc (M Int)
    offsets = accM

    mat :: Acc (Array DIM3 Double)
    mat = generate (lift (Z :. M.nrows m :. M.ncols m :. sz)) $
      \ ix ->
        let (r, c, _) = unlift $ unindex3 ix :: (Exp Int, Exp Int, Exp Int)
            ix'       = index2 r c
        in
        cond (accM ! ix' == 0)
             0
             liftedZ
    liftedZ = lift z

matrixMatrixlet :: forall a b . (Elt b, Num a) => Exp b -> Acc (Matrixlet a) -> (Exp DIM2 -> Exp b) -> Acc (Matrixlet b)
matrixMatrixlet zero t k =
  lift (sz, indices, offsets, generate (shape mat) k')
  where
    (sz, indices, offsets, mat) = unlift t
                                    :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    sz' = the sz

    k' :: Exp DIM3 -> Exp b
    k' ix =
      cond (indices ! ix')
           (arrayArraylet sz off arrayletIx $ \ ix ->
              let (r', c') = unlift $ unindex2 ix
              in k $ index2 (r * sz' + r') (c * sz' + c'))
           zero
      where
        ix'                = index2 r c :: Exp DIM2
        (r, c, arrayletIx) = unlift $ unindex3 ix :: (Exp Int, Exp Int, Exp Int)
        off                = offsets ! ix'

-- Assumes there is always one value on every column
foldColsMatrixlet :: forall a . (Elt a) => (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldColsMatrixlet f mat = fold1 f $ fold1 f arr
  where
    (_, _, _, arr) = unlift mat
                       :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))


-- | Runs Accelerate computation on GPU
fromAcc :: Acc (V Bool) -> U.Vector Bool
fromAcc = U.map word8ToBool . U.convert . toVectors . run
  where
    word8ToBool 0 = False
    word8ToBool 1 = True
    word8ToBool n = error $ "word8ToBool: Got binary argument '" P.++ show n P.++ "'"

toAcc :: DIM1 -> U.Vector Double -> Acc (V Double)
toAcc sh = lift . fromVectors sh . U.convert

hard' :: Exp Double -> Exp Bool
hard' = (> 0)

ldpc :: Acc (Matrixlet Double) -> Int -> Acc (V Double) -> Acc (V Bool)
ldpc mLet maxIterations orig_lam = {- traceShow msg $ -} map hard' $ loop 0 mLet orig_lam
  where
    loop :: Int -> Acc (Matrixlet Double) -> Acc (V Double) -> Acc (V Double)
    loop !n ne lam =
      let (finalN, _, r) = unlift (awhile loopCond loopBody liftedInit)
                            :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
      in
      acond (the finalN >= lift maxIterations)
            (lift orig_lam)
            (lift lam)
      where
        liftedInit :: Acc (Scalar Int, Matrixlet Double, V Double)
        liftedInit = lift (unit (lift n), ne, lam)

    loopCond :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Bool)
    loopCond t = unit $ not (the (all (== lift False) ans) || n >= lift maxIterations)
      where
        (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
        n              = the n0

        ans :: Acc (V Bool)
        ans = foldColsMatrixlet (/=) $ matrixMatrixlet (lift False) mLet
          $ \ ix ->
                let (r, c) = unlift $ unindex2 ix :: (Exp Int, Exp Int)
                in hard' (lam ! lift (Z :. c))

    loopBody :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Int, Matrixlet Double, V Double)
    loopBody  = undefined
      -- where
      --   ne_tanh'mat :: Matrixlet Double
      --   ne_tanh'mat =
      --     imapMatrixlet
      --       (\ (m,n) v -> tanh (- ((lam U.! n - v) / 2)))
      --       ne

      --   ne_tanhMulted :: V StableDiv
      --   ne_tanhMulted = foldColsMatrixletU (\_ v -> lit v) smult ne_tanh'mat

      --   ne' :: Matrixlet Double
      --   ne' =
      --     imapMatrixlet
      --       (\ (m,n) v -> -2 * atanh' ((ne_tanhMulted U.! m) `sdiv` v))
      --       ne_tanh'mat

      --   lam' :: V Double
      --   lam' = U.zipWith (+) orig_lam $ foldRowsMatrixlet (+) ne'

