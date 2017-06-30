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

type Matrixlet2D a = Array DIM2 a

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
foldColsMatrixlet f = foldRowsMatrixlet f . transposeMatrixlet
-- foldColsMatrixlet f mat =
--   generate (lift (Z :. mletRowCount*the sz)) $ \ix ->
--     let r = unindex1 $ unlift ix
--         (rd, -- Block row
--          rr  -- Row inside block
--          )    = r `divMod` the sz
--     in
--     A.snd $
--     while ((< mletColCount) . A.fst)
--           (\p ->
--             let (c, v) = unlift p :: (Exp Int, Exp a)
--                 x      = m ! lift (Z :. rd :. c :. rr)
--             in
--             lift (c+1, f x v)
--           )
--           (lift (1::Int, m ! lift (Z :. rd :. (0 :: Int) :. rr)))
--   where
--     (sz, indices, offsets, m)  = unlift mat
--         :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
--     (mletRowCount, mletColCount, arrayletSize) = unlift $ unindex3 (shape m) :: (Exp Int, Exp Int, Exp Int)

foldRowsMatrixlet :: forall a . Elt a => (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldRowsMatrixlet f mat =
  generate (lift (Z :. mletColCount*the sz)) $ \ix ->
    let c = unindex1 $ unlift ix
        (cd, -- Block column
         cr  -- Column inside block
         )    = c `divMod` the sz
        arrRow = (cr + (the sz `div` 2)) `mod` the sz
    in
    sfoldl f
           (m ! (lift (Z :. cd :. arrRow :. (0 :: Int))))
           (lift (Z :. cd :. arrRow))
           (A.tail swappedM)

    -- A.snd $
    -- while ((< mletRowCount) . A.fst)
    --       (\p ->
    --         let (r, v) = unlift p :: (Exp Int, Exp a)
    --             x      = m ! lift (Z :. r :. cd :. cr)
    --         in
    --         lift (r+1, f x v)
    --       )
    --       (lift (1::Int, m ! lift (Z :. (0 :: Int) :. cd :. cr)))
  where
    swappedM = backpermute
      (lift (Z :. mletColCount :. arrayletSize :. mletRowCount))
      (\p ->
        let (r, c, i) = unlift $ unindex3 p :: (Exp Int, Exp Int, Exp Int)
        in
        lift (Z :. c :. i :. r))
      m
    -- (r, c, aLetIx) --> (r, aLetIx, c)
    -- swappedM = backpermute
    --   (lift (Z :. mletRowCount :. arrayletSize :. mletColCount))
    --   (\p ->
    --     let (r, c, i) = unlift $ unindex3 p :: (Exp Int, Exp Int, Exp Int)
    --     in
    --     lift (Z :. r :. i :. c))
    --   m

    (sz, indices, offsets, m) = unlift mat --unlift $ transposeMatrixlet mat
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    (mletRowCount, mletColCount, arrayletSize) = unlift $ unindex3 (shape m) :: (Exp Int, Exp Int, Exp Int)
    colCount = arrayletSize * mletColCount

-- -- TODO: Check indexing transformation
-- dim3ToDim2 :: (P.Num a, P.Integral a) => a -> a -> a -> a -> a -> (a, a)
-- dim3ToDim2 sz offset r c arrayletIx =
--     ((sz*r) + arrayletIx, (sz*c) + ((arrayletIx + offset) `mod` sz))

imapMatrixlet :: forall a b . (Elt a, Elt b) => ((Exp Int, Exp Int) -> Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
imapMatrixlet f t = lift (sz, indices, offsets, imap f' mat)
  where
    arrayletIxTransform :: Exp Int -> Exp Int -> ((Exp Int, Exp Int) -> r) -> r
    arrayletIxTransform offset arrayletIx g =
      g (arrayletIx, (arrayletIx + offset) `mod` the sz)

    f' :: Exp DIM3 -> Exp a -> Exp b
    f' ix x =
      arrayletIxTransform (offsets ! ix') arrayletIx $ \ (r', c') ->
        let r'' = (r*the sz) + r'
            c'' = (c*the sz) + c'
        in
        f (r'', c'') x
      where
        (r, c, arrayletIx) = unlift $ unindex3 ix :: (Exp Int, Exp Int, Exp Int)
        ix' = index2 r c :: Exp DIM2

-- imapMatrixlet f t = lift (sz, indices, offsets, imap f' mat)
--   where
--     f' :: Exp DIM3 -> Exp a -> Exp b
--     f' ix = f (dim3ToDim2 (the sz) (offsets ! ix') r c arrayletIx)
--       where
--         (r, c, arrayletIx) = unlift $ unindex3 ix :: (Exp Int, Exp Int, Exp Int)
--         ix' = index2 r c :: Exp DIM2

    (sz, indices, offsets, mat) = unlift t
                                    :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

mapMatrixlet ::  forall a b . (Elt a, Elt b) => (Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
mapMatrixlet f t = lift (sz, indices, offsets, map f mat)
  where
    (sz, indices, offsets, mat) = unlift t
                                    :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))

foldMapColsMatrixlet
          :: (Elt a, Elt b)
             => (Exp a -> Exp b)
             -> (Exp b -> Exp b -> Exp b)
             -> Acc (Matrixlet a)
             -> Acc (V b)
foldMapColsMatrixlet f g t = foldColsMatrixlet' g $ mapMatrixlet f t

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

foldColsMatrixlet' :: forall a . (Elt a) => (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldColsMatrixlet' f mLet = fold1 f m2D
  where
    m2D = to2DArray mLet

foldRowsMatrixlet' :: forall a . (Elt a) => (Exp a -> Exp a -> Exp a) -> Acc (Matrixlet a) -> Acc (V a)
foldRowsMatrixlet' f mLet = fold1 f $ transpose m2D
  where
    (sz, indices, offsets, m) = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    m2D = to2DArray mLet

foldRowsMatrixlet'' :: forall a . (Elt a) => (Exp a -> Exp a -> Exp a) -> Exp a -> Acc (Matrixlet a) -> Acc (V a)
foldRowsMatrixlet'' f z mLet =
  permute f (fill (lift (Z :. colCount)) z) reindex m
  where
    reindex :: Exp DIM3 -> Exp DIM1
    reindex p = index1 ((i*the sz) + k)
      where
        (i, j, k) = unlift $ unindex3 p :: (Exp Int, Exp Int, Exp Int)

    (sz, indices, offsets, m) = unlift (transposeMatrixlet mLet)
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    (colCount, rowCount, _) = unlift . unindex3 $ shape m :: (Exp Int, Exp Int, Exp Int)

imapMatrixlet' :: forall a b . (Elt a, Elt b) => ((Exp Int, Exp Int) -> Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
imapMatrixlet' f mLet =
  overMatrixlet mLet (\m ->
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
to2DArray :: forall a . Elt a => Acc (Matrixlet a) -> Acc (Array DIM2 a)
to2DArray mLet =
  backpermute shape2D (from2D (the sz)) m
  where
    (sz, indices, offsets, m) = unlift mLet
        :: (Acc (Scalar Int), Acc (M Bool), Acc (M Int), Acc (Array DIM3 a))
    (mLetRowCount, mLetColCount, _) = unlift $ unindex3 $ shape m
      :: (Exp Int, Exp Int, Exp Int)
    shape2D =
      lift (Z :. mLetRowCount*the sz :. mLetColCount*the sz)
        :: Exp DIM2
-- | Abstracts of internal 3D index representation
overMatrixlet :: forall a b . (Elt a, Elt b) =>
  Acc (Matrixlet a) -> (Acc (Matrixlet2D (DIM2, a)) -> Acc (Matrixlet2D b)) -> Acc (Matrixlet b)
overMatrixlet mLet f =
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
          :: Acc (Matrixlet2D (DIM2, a)))

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

        ans :: Acc (V Bool)
        ans =
          A.map A.even $
          foldColsMatrixlet' (+) $
          imapMatrixlet' (\ (r, c) v ->
            cond (c_hat !! c && v /= 0)
                 1
                 0 :: Exp Int) $
          mLet

          -- foldColsMatrixlet (/=) $
          -- imapMatrixlet' (\ (r, c) _ -> c_hat !! r) $
          -- mLet

        c_hat :: Acc (V Bool)
        c_hat = map hard' lam

    loopBody :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Int, Matrixlet Double, V Double)
    loopBody t = lift (unit (n + 1) :: Acc (Scalar Int), ne', lam')
      where
        (!n0, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))

        n = the n0 :: Exp Int

        ne_tanh'mat :: Acc (Matrixlet Double)
        ne_tanh'mat =
          imapMatrixlet'
            (\ (m,n) v -> tanh (- ((lam ! (lift (Z :. n)) - v) / 2)))
            ne

        ne_tanhMulted :: Acc (V StableDiv)
        ne_tanhMulted = foldMapColsMatrixlet lit smult ne_tanh'mat

        ne' :: Acc (Matrixlet Double)
        ne' =
          imapMatrixlet'
            (\ (m,n) v -> -2 * atanh'' (ne_tanhMulted ! (lift (Z :. m)) `sdiv` v))
            ne_tanh'mat

        lam' :: Acc (V Double)
        lam' = zipWith (+) orig_lam $ foldRowsMatrixlet'' (+) 0 ne'

