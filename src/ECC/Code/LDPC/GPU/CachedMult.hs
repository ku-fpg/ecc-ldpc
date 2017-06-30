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

type Arraylet a =
  (,)
    (Scalar Int)
    (V a)

mkArraylet :: Elt a => Acc (Scalar Int) -> Acc (V a) -> Acc (Arraylet a)
mkArraylet sz v = lift (sz, v)

arrayArraylet :: Acc (Scalar Int) -> Exp Int -> Exp Int -> (Exp DIM2 -> Exp b) -> Exp b
arrayArraylet sz0 off arrayletIx k =
  k (index2 arrayletIx ((arrayletIx + off) `mod` sz))
  where
    sz = the sz0

foldRowsArraylet :: forall a . (Elt a) => Acc (Arraylet a) -> Acc (V a)
foldRowsArraylet p = b ++ a
  where
    (sz0, m) = unlift p :: (Acc (Scalar Int), Acc (V a))
    sz = the sz0
    splitDist = (A.length m - sz)
    a = take splitDist m
    b = drop splitDist m

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


transposeMatrixlet :: forall a . (Elt a) => Acc (Matrixlet a) -> Acc (Matrixlet a)
transposeMatrixlet p =
  lift
    (sz
    ,transpose indices
    ,transpose offsets
    ,backpermute (lift (Z :. origC :. origR :. arrayletSize)) go m
    )
  where
    (origR, origC, arrayletSize) =
      unlift $ unindex3 $ shape m :: (Exp Int, Exp Int, Exp Int)

    go q =
      let (r, c, ix) = unlift $ unindex3 q :: (Exp Int, Exp Int, Exp Int)
          -- This should be equivalent to swapping the two halves of the
          -- arraylet:
          ix'        = (ix + (arrayletSize `quot` 2)) `mod` arrayletSize
      in
      lift (Z :. c :. r :. ix')

    (sz, indices, offsets, m) = unlift p
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

foldColsMatrixlet :: forall a . (Elt a, P.Num a, Lift Exp a, a ~ Plain a) =>
  (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldColsMatrixlet f mLet = fold1 f m2D
  where
    m2D = to2DArray mLet

foldColsMatrixlet' :: forall a . (Elt a) =>
  Exp a -> (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldColsMatrixlet' z f mLet = fold1 f m2D
  where
    m2D = to2DArrayWith z mLet

foldRowsMatrixlet :: forall a . (Elt a, P.Num a, Lift Exp a, a ~ Plain a) =>
  (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldRowsMatrixlet f mLet = fold1 f $ transpose m2D
  where
    (sz, indices, offsets, m) = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    m2D = to2DArray mLet

imapMatrixlet :: forall a b . (Elt a, P.Num a, a ~ Plain a, Lift Exp a, Elt b) =>
  ((Exp Int, Exp Int) -> Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
imapMatrixlet f mLet =
  matrixletAs2D
    mLet
    (\m ->
      lift (imap (\p -> f' p . A.snd) m))
  where
    f' :: Exp DIM2 -> Exp a -> Exp b
    f' p x =
      let (r, c) = unlift $ unindex2 p :: (Exp Int, Exp Int)
      in
      f (r, c) x

    m2D = to2DArray mLet
    (sz, indices, offsets, _) = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

mapMatrixlet :: forall a b . (Elt a, Elt b) => (Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
mapMatrixlet f = overMatrixlet (A.map f)


-- TODO: Verify that these transformations are correct.
from2D :: Exp Int -> Exp DIM2 -> Exp DIM3
from2D sz p =
  let (r, c) = unlift $ unindex2 p :: (Exp Int, Exp Int)
  in
  lift
    (Z
    :. (r `div` sz)
    :. (c `div` sz)
    :. (r `mod` sz)
    )

to2D :: Exp Int -> Exp DIM3 -> Exp DIM2
to2D sz p =
  let (r, c, i) = unlift $ unindex3 p :: (Exp Int, Exp Int, Exp Int)
      i'        = (i + (sz `quot` 2)) `mod` sz
  in
  lift
    (Z
    :. ((r*sz)+i)
    :. ((c*sz)+i')
    )

to2DArray :: forall a . (Lift Exp a, Elt a, P.Num a, a ~ Plain a) =>
  Acc (Matrixlet a) -> Acc (Array DIM2 a)
to2DArray = to2DArrayWith (lift (0 :: a))

to2DArrayWith :: forall a . (Elt a) =>
  Exp a -> Acc (Matrixlet a) -> Acc (Array DIM2 a)
to2DArrayWith z mLet =
  generate shape2D $ \rc ->
    let (r, c)                    = unlift $ unindex2 rc :: (Exp Int, Exp Int)
        mLetRCI                   = from2D (the sz) rc
        (mLetR, mLetC, mLetArrIx) = unlift $ unindex3 mLetRCI
                                      :: (Exp Int, Exp Int, Exp Int)
        mLetRC                    = lift (Z :. mLetR :. mLetC) :: Exp DIM2
    in
    cond (indices ! mLetRC)
         (m ! mLetRCI)
         z
  where
    nonzeroMLetIndices :: Acc (V DIM2)
    nonzeroMLetIndices =
      A.map A.fst $
      flatten $
      indexed indices

    (sz, indices, offsets, m) = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    (mLetRowCount, mLetColCount, _) = unlift $ unindex3 $ shape m
      :: (Exp Int, Exp Int, Exp Int)
    shape2D =
      lift (Z :. mLetRowCount*the sz :. mLetColCount*the sz)
        :: Exp DIM2
-- | Abstracts of internal 3D index representation
matrixletAs2D :: forall a b . (Elt a, Elt b) =>
  Acc (Matrixlet a) -> (Acc (M (DIM2, a)) -> Acc (M b)) -> Acc (Matrixlet b)
matrixletAs2D mLet f =
    lift
      (sz
      ,indices
      ,offsets
      ,backpermute (shape m) (to2D (the sz)) m'
      )
  where
    (sz, indices, offsets, m)  = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

    (mLetRowCount, mLetColCount, _) = unlift $ unindex3 $ shape m
      :: (Exp Int, Exp Int, Exp Int)

    dropDim :: Exp DIM3 -> Exp DIM2
    dropDim p =
      let (r, c, _) = unlift $ unindex3 p :: (Exp Int, Exp Int, Exp Int)
      in
      lift (Z :. r :. c)

    dropDim' :: Exp (DIM3, a) -> Exp (DIM2, a)
    dropDim' p =
      let (d, v) = unlift p :: (Exp DIM3, Exp a)
      in
      lift (dropDim d, v)


    shape2D =
      lift (Z :. mLetRowCount*the sz :. mLetColCount*the sz)
        :: Exp DIM2

    m' =
      f (backpermute shape2D (from2D (the sz)) (A.map dropDim' (indexed m))
          :: Acc (M (DIM2, a)))

overMatrixlet :: forall a b . (Elt a, Elt b) =>
   (Acc (Array DIM3 a) -> Acc (Array DIM3 b)) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
overMatrixlet f mLet = lift (sz, indices, offsets, f m)
  where
    (sz, indices, offsets, m)  = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

-- | Do the given coordinates give a non-zero value in a circulant with the
-- given size and offset parameters? The coordinates are given in arraylet
-- coodinate space.
nonzeroArrayletIx :: Exp Int -> Exp Int -> Exp DIM2 -> Exp Bool
nonzeroArrayletIx sz offset d =
  let (r, c) = unlift $ unindex2 d :: (Exp Int, Exp Int)
  in
  c == (r + offset) `mod` sz

ldpc :: Acc (Matrixlet Double) -> Int -> Acc (V Double) -> Acc (V Bool)
ldpc mLet maxIterations orig_lam = {- traceShow msg $ -} map hard' $ loop 0 mLet orig_lam
  where
    loop :: Int -> Acc (Matrixlet Double) -> Acc (V Double) -> Acc (V Double)
    loop !n ne lam =
      let (finalN, ne, r) = unlift (awhile loopCond loopBody liftedInit)
                            :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
      in
      acond (the finalN >= lift maxIterations)
            (lift orig_lam) --(lift r)--(lift orig_lam)
            (lift r)
      where
        liftedInit :: Acc (Scalar Int, Matrixlet Double, V Double)
        liftedInit = lift (unit (lift n), ne, lam)

    loopCond :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Bool)
    loopCond t = unit $ (n < lift maxIterations) && not (the (A.any (== lift False) ans))
      where
        (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
        n              = the n0
        (sz, indices, _, m)  = unlift ne
            :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 Double))
        (mLetRows, mLetCols, _) = unlift $ unindex3 $ shape m
          :: (Exp Int, Exp Int, Exp Int)

        ans :: Acc (V Bool)
        ans =
          fold1 (/=) $
          A.map (\p ->
            let (d, v) = unlift p :: (Exp DIM3, Exp Double)
                (_, c) = unlift . unindex2 $ to2D (the sz) d
                  :: (Exp Int, Exp Int)
            in
            c_hat !! c) $
          generate (lift (Z :. mLetRows*the sz :. mLetCols*the sz)) (\ix ->
            undefined)

        w :: Acc (V (DIM3, Double))
        w =
          A.afst $
          A.filter (\p ->
            let (ix, v) = unlift p :: (Exp DIM3, Exp Double)
                (r, c, i) = unlift $ unindex3 ix :: (Exp Int, Exp Int, Exp Int)
            in
            indices ! lift (Z :. r :. c)) $
          indexed $
          m

          -- foldColsMatrixlet (/=) $
          -- imapMatrixlet (\ (r, c) _ -> c_hat !! c) $
          -- mLet

          -- map A.even $
          -- permute
          --   (+)
          --   (fill (shape c_hat) 0)
          --   (\mLetRC ->
          --     let (mLetR, mLetC, _) = unlift $ unindex3 mLetRC
          --           :: (Exp Int, Exp Int, Exp Int)
          --         (r, c) = to2D mLetRC
          --     in
          --     c
          --   )
          --   (map (\b ->
          --     cond b (lift (1 :: Int)) (lift (0 :: Int)))
          --         m')

          -- generate (lift (shape c_hat)) $ \ix ->
          --   lift False

        (_, _, _, m') = unlift $ imapMatrixlet (\ (r, c) _ -> c_hat !! c) mLet
            :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 Bool))

        c_hat :: Acc (V Bool)
        c_hat = map hard' lam

    loopBody :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Int, Matrixlet Double, V Double)
    loopBody t = lift (unit (n + 1) :: Acc (Scalar Int), ne', lam')
      where
        (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))

        n = the n0 :: Exp Int

        ne_tanh'mat :: Acc (Matrixlet Double)
        ne_tanh'mat =
          imapMatrixlet
            (\ (m,n) v -> tanh (- ((lam ! (lift (Z :. n)) - v) / 2)))
            ne

        ne_tanhMulted :: Acc (V StableDiv)
        ne_tanhMulted =
          foldColsMatrixlet' (lit 1) smult $
          imapMatrixlet (const lit) ne_tanh'mat

        ne' :: Acc (Matrixlet Double)
        ne' =
          imapMatrixlet
            (\ (m,n) v ->
              -2 * atanh'' ((ne_tanhMulted ! (lift (Z :. m))) `sdiv` v))
            ne_tanh'mat

        lam' :: Acc (V Double)
        lam' = zipWith (+) orig_lam $ foldRowsMatrixlet (+) ne'

