{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, DeriveFunctor #-}
{-# LANGUAGE FlexibleContexts, MultiWayIf, TypeFamilies #-}
module ECC.Code.LDPC.GPU.CachedMult where

-- Uses the StableDiv data structure to cache multiplications.
-- Uses Accelerate to run on the GPU

import Prelude hiding ((==), (/=), (>=), (<), (>), all, map, (||), (&&), not, Num, snd, zipWith, (++), length, take, drop, RealFloat, Eq, Fractional, Floating, (!!))
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

import Data.Foldable (foldl')

import Data.Array.Accelerate as A
import Data.Array.Accelerate.IO
import Data.Array.Accelerate.LLVM.PTX

import Control.Monad.State hiding (lift)
import Data.Traversable (for)
import Data.Monoid

code :: Code
code = mkLDPC_Code "gpu-arraylet-cm" E.encoder decoder

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a rate maxIterations orig_lam =
  Just $ fromAcc (ldpc mLet maxIterations (toAcc sh orig_lam))
  where mLet = initMatrixlet 0 a
        sh   = Z :. U.length orig_lam

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


{-# INLINE atanh'' #-}
atanh'' :: (RealFloat a, Eq a, Fractional a, Floating a) => Exp a -> Exp a
atanh'' x =
  cond (x == 1 || x == -1)
       (signum x * 18.714973875118524)
       (atanh x)

---------------------------------------------------------------------


type StableDiv
  = ((,)
       Double  -- | The "worst" value for division: the closest to zero
       Double) -- | The result of the multiplications,
                   --   excluding the "worst" value

absMinMax :: Exp Double -> Exp Double -> Exp (Double, Double)
absMinMax x y =
  cond (abs x < abs y)
       (lift (x, y))
       (lift (y, x))

lit :: Exp Double -> Exp StableDiv
lit x =
  cond (x >= 1)
       (lift (1 :: Double, x))
       (lift (x, 1 :: Double))

smult :: Exp StableDiv -> Exp StableDiv -> Exp StableDiv
smult p q =
    lift (minOfMins, (b * maxOfMins * d))
    where
      minOfMins, maxOfMins :: Exp Double
      (minOfMins, maxOfMins) = unlift (absMinMax a c)

      (a, b) = unlift p
      (c, d) = unlift q

sdiv :: Exp StableDiv -> Exp Double -> Exp Double
sdiv p c =
  cond (a == c)
       b
       (a * (b/c))
  where
    (a, b) = unlift p

type V a = Array DIM1 a
type M a = Array DIM2 a


type Arraylet  a =
  (,)
    Int   -- | Offset
    (V a)

type Matrixlet a =
  (,,,)
    (Scalar Int) -- | Arraylet size
    (M Int)      -- | Arraylet indices (given as column numbers)
    (V Int)      -- | Arraylet offsets
    (M a)        -- | Arraylets (one arraylet per column)



initMatrixlet :: forall a. (U.Unbox a, Elt a, a ~ Plain a, Lift Exp a) =>
  a -> Q.QuasiCyclic Integer -> Acc (Matrixlet a)
initMatrixlet zero (Q.QuasiCyclic sz qm) =
    lift (unit (lift sz), accIndices, offsets, arraylets)
  where
    indices :: M.Matrix Int
    indices =
      flip evalState 0 $
      for qm $ \n ->
        case n of
          0 -> pure 0
          _ -> modify (+1) *> get

    accIndices :: Acc (M Int)
    accIndices =
      use $
      fromVectors (Z :. M.nrows qm :. M.ncols qm) $
      U.convert $
      M.getMatrixAsVector indices

    -- Arraylet matrix (initially all 'zero's)
    arraylets :: Acc (M a)
    arraylets =
      fill (lift (Z
                  :. sz
                  :. (M.nrows qm*M.ncols qm)-zeroArrayletCount))
           (lift zero)

    -- Offset vector
    offsets :: Acc (V Int)
    offsets =
      use $
      fromList (Z :. (M.nrows qm*M.ncols qm)-zeroArrayletCount) $
      foldMap (\n ->
        case n of
          0 -> []
          _ -> [g n]) $
      qm

    zeroArrayletCount :: Int
    zeroArrayletCount =
      getSum $
      foldMap (\n ->
        case n of
          0 -> 1
          _ -> 0)
      qm

    -- The number must be a power of two, because there is only one bit set.
    g :: Integer -> Int
    g x | x `testBit` 0 = 0
        | x P.== 0      = error "got to zero; should never happen"
        | otherwise     = 1 + g (x `shiftR` 1)

-- | Do the given coordinates give a non-zero value in a circulant with the
-- given size and offset parameters? The coordinates are given in arraylet
-- coodinate space.
nonzeroArrayletIx :: Exp Int -> Exp Int -> Exp DIM2 -> Exp Bool
nonzeroArrayletIx sz offset d =
  let (r, c) = unlift $ unindex2 d :: (Exp Int, Exp Int)
  in
  c == (r + offset) `mod` sz

ldpc :: Acc (Matrixlet Double) -> Int -> Acc (V Double) -> Acc (V Bool)
ldpc = error "Not implemented"




-- ldpc :: Acc (Matrixlet Double) -> Int -> Acc (V Double) -> Acc (V Bool)
-- ldpc mLet maxIterations orig_lam = {- traceShow msg $ -} map hard' $ loop 0 mLet orig_lam
--   where
--     loop :: Int -> Acc (Matrixlet Double) -> Acc (V Double) -> Acc (V Double)
--     loop !n ne lam =
--       let (finalN, ne, r) = unlift (awhile loopCond loopBody liftedInit)
--                             :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
--       in
--       acond (the finalN >= lift maxIterations)
--             (lift orig_lam) --(lift r)--(lift orig_lam)
--             (lift r)
--       where
--         liftedInit :: Acc (Scalar Int, Matrixlet Double, V Double)
--         liftedInit = lift (unit (lift n), ne, lam)

--     loopCond :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Bool)
--     loopCond t = unit $ (n < lift maxIterations) && not (the (A.any (== lift False) ans))
--       where
--         (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
--         n              = the n0
--         (sz, indices, _, m)  = unlift ne
--             :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 Double))
--         (mLetRows, mLetCols, _) = unlift $ unindex3 $ shape m
--           :: (Exp Int, Exp Int, Exp Int)

--         ans :: Acc (V Bool)
--         ans = undefined

--         c_hat :: Acc (V Bool)
--         c_hat = map hard' lam

--     loopBody :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Int, Matrixlet Double, V Double)
--     loopBody t = lift (unit (n + 1) :: Acc (Scalar Int), ne', lam')
--       where
--         (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))

--         n = the n0 :: Exp Int

--         ne_tanh'mat :: Acc (Matrixlet Double)
--         ne_tanh'mat =
--           imapMatrixlet
--             (\ (m,n) v -> tanh (- ((lam ! (lift (Z :. n)) - v) / 2)))
--             ne

--         ne_tanhMulted :: Acc (V StableDiv)
--         ne_tanhMulted =
--           foldColsMatrixlet' (lit 1) smult $
--           imapMatrixlet (const lit) ne_tanh'mat

--         ne' :: Acc (Matrixlet Double)
--         ne' =
--           imapMatrixlet
--             (\ (m,n) v ->
--               -2 * atanh'' ((ne_tanhMulted ! (lift (Z :. m))) `sdiv` v))
--             ne_tanh'mat

--         lam' :: Acc (V Double)
--         lam' = zipWith (+) orig_lam $ foldRowsMatrixlet (+) ne'

