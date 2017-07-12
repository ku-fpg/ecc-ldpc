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
import qualified Data.Vector as V
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as S
import qualified Data.Matrix.QuasiCyclic as Q
import Debug.Trace

import qualified ECC.Code.LDPC.Fast.Encoder as E

import Data.Foldable (foldl')

import Data.Array.Accelerate as A
import Data.Array.Accelerate.IO
import Data.Array.Accelerate.LLVM.PTX

import Control.Monad.State hiding (lift)
import Data.Traversable (for)
import Data.Monoid hiding (All, Any)

code :: Code
code = ldpc `seq` mkLDPC_Code "gpu-arraylet-cm" E.encoder decoder

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder arr@(Q.QuasiCyclic sz qm) = a `seq` b `seq` c `seq` d `seq` e `seq` f `seq` \rate maxIterations orig_lam ->
  let sh   = Z :. U.length orig_lam
  in
  Just $ U.map word8ToBool $ U.convert $ toVectors $ ldpc (a, b, c, d, e, f, run $ unit $ lift maxIterations, fromVectors sh (U.convert orig_lam))
  where
    (a, b, c, d, e, f) = run $ initMatrixlet (unit $ lift (0 :: Double)) (unit (lift sz)) (unit $ lift zeroArrayletCount) (use nonzeroBlocks) (unit $ lift qmSh) (use arrayletIndices) (use offsets)
    qmSh = Z :. M.nrows qm :. M.ncols qm

    -- The number must be a power of two, because there is only one bit set.
    g :: Integer -> Int
    g x | x `testBit` 0 = 0
        | x P.== 0      = error "got to zero; should never happen"
        | otherwise     = 1 + g (x `shiftR` 1)

    -- Offset vector
    offsets :: V Int
    offsets =
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
          _ -> 0) $
      qm

      -- How many zero entries to skip
    arrayletIndices :: M Int
    arrayletIndices =
      fromVectors qmSh $
      U.convert $
      (flip evalState 0
        . V.mapM (\x ->
            if x P.== 0
            then do
              y <- get
              modify (+1)
              return y
            else get)
          ) $
      M.getMatrixAsVector qm


    nonzeroBlocks :: M Bool
    nonzeroBlocks =
      fromVectors qmSh
        (S.map boolToWord8
          (U.convert (M.getMatrixAsVector (fmap (P./= 0) qm))))



-- | Compiled ldpc'
ldpc :: (Scalar Int, Scalar DIM2, V DIM2, V Int, M Int, M Double, Scalar Int, V Double) -> V Bool
ldpc =
  -- traceShow computation $
  run1 (\t ->
    let (a, b, c, d, e, f, maxIterations, orig_lam) = unlift t ::
          (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M Double), Acc (Scalar Int), Acc (V Double))
    in
    ldpc' (lift (a, b, c, d, e, f)) maxIterations orig_lam)

-- | Runs Accelerate computation on GPU
fromAcc :: Acc (V Bool) -> U.Vector Bool
fromAcc = U.map word8ToBool . U.convert . toVectors . run

word8ToBool :: Word8 -> Bool
word8ToBool 0 = False
word8ToBool 1 = True
word8ToBool n = error $ "word8ToBool: Got binary argument '" P.++ show n P.++ "'"

boolToWord8 :: Bool -> Word8
boolToWord8 False = 0
boolToWord8 True  = 1

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
  (,,,,,)
    (Scalar Int)  -- | Arraylet size
    (Scalar DIM2) -- | "Actual" matrix size
    (V DIM2)      -- | Arraylet block coordinates
    (V Int)       -- | Arraylet offsets
    (M Int)       -- | Arraylet indices
    (M a)         -- | Arraylets (one arraylet per column)

initMatrixlet :: forall a. (U.Unbox a, Elt a, a ~ Plain a, Lift Exp a) =>
  Acc (Scalar a) -> Acc (Scalar Int) -> Acc (Scalar Int) -> Acc (M Bool) -> Acc (Scalar DIM2) -> Acc (M Int) -> Acc (V Int) -> Acc (Matrixlet a)
initMatrixlet zero sz zeroArrayletCount nonzeroBlocks sh arrayletIndices offsets =
    lift (sz, dim, blockCoords, offsets, arrayletIndices, arraylets)
  where
    (rowCount, colCount) = unlift $ unindex2 $ the sh :: (Exp Int, Exp Int)
    dim = unit $ lift (Z :. rowCount*the sz :. colCount*the sz)

    blockCoords :: Acc (V DIM2)
    blockCoords =
      A.map A.fst $
      A.afst $
      A.filter A.snd $
      indexed $
      nonzeroBlocks

    -- Arraylet matrix (initially all 'zero's)
    arraylets :: Acc (M a)
    arraylets =
      fill (lift (Z
                  :. the sz
                  :. (rowCount*colCount)-the zeroArrayletCount))
           (lift (the zero))

-- | Do the given coordinates give a non-zero value in a circulant with the
-- given size and offset parameters? The coordinates are given in arraylet
-- coodinate space.
nonzeroArrayletIx :: Exp Int -> Exp Int -> Exp DIM2 -> Exp Bool
nonzeroArrayletIx sz offset d =
  let (r, c) = unlift $ unindex2 d :: (Exp Int, Exp Int)
  in
  c == (r + offset) `mod` sz

foldRowsMatrixlet :: forall a . (Elt a, a ~ Plain a) => (Exp a -> Exp a -> Exp a) -> Exp a -> Acc (Matrixlet a) -> Acc (V a)
foldRowsMatrixlet f z mLet =
  permute
    (\p q ->
      let a = unlift p :: Exp a
          b = unlift q :: Exp a
      in
      f a b)
    (A.replicate (lift (Z :. realColCount)) (unit z))
    (\ix ->
      let (r0, c0)         = unlift $ unindex2 ix :: (Exp Int, Exp Int)
          (blockR, blockC) = unlift $ unindex2 $ blockCoords !! c0 :: (Exp Int, Exp Int)
          offset           = offsets !! c0
      in index1 ((blockC*the sz) + ((r0 + offset) `mod` the sz))) $
    arraylets
  where
    (sz, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))
    arrayletsTr = transpose arraylets

    (_, realColCount) = unlift $ unindex2 $ the realDim :: (Exp Int, Exp Int)

foldColsMatrixlet :: forall a . (Elt a, a ~ Plain a) => (Exp a -> Exp a -> Exp a) -> Exp a -> Acc (Matrixlet a) -> Acc (V a)
foldColsMatrixlet f z mLet =

  -- generate (lift (Z :. realRowCount)) $ \ix0 ->
  --   let r' = unindex1 ix0
  --       r  = r' `mod` the sz
  --   in
  --   sfoldl
  --     f
  --     z
  --     (lift Z)

  -- A.asnd $
  -- awhile (unit . (< realRowCount) . the . A.afst) (\p ->
  --   let (r'0, v) = unlift p :: (Acc (Scalar Int), Acc (Scalar a))
  --       r' = the r'0
  --       r  = r' `mod` the sz
  --   in
  --   lift
  --   (unit $ r'+1
  --   ,unit $ fold undefined undefined
  --       (slit
  --         (r' `div` the sz)
  --         (the sz)
  --         (slice arraylets (lift (Z :. r :. All))))
  --   ))
  --   (lift (0, z))

  fold f z $
  permute
    (\p q ->
      let a = unlift p :: Exp a
          b = unlift q :: Exp a
      in
      a)
    (A.replicate (lift (Z :. realRowCount :. blockColCount)) (unit z))
    (\ix ->
      let (r0, c0)         = unlift $ unindex2 ix :: (Exp Int, Exp Int)
          (blockR, blockC) = unlift $ unindex2 $ blockCoords !! c0 :: (Exp Int, Exp Int)
          offset           = offsets !! c0
      in lift ((Z :. (blockR*the sz) + r0) :. blockC)) $
    arraylets

  -- fold f z $
  -- -- transpose $
  -- backpermute
  --   (lift (Z :. realRowCount :. blockColCount))
  --   (\ix ->
  --     let (r', c') = unlift $ unindex2 ix :: (Exp Int, Exp Int)
  --         -- arrayIx  = arrayletIndices ! lift (Z :. (r' `div` the sz) :. c')
  --         r        = r' `mod` the sz
  --         c        = c' + ((r' `div` the sz) * blockColCount)
  --     in
  --     (lift (Z :. r :. c))
  --   )
  --   arraylets

    -- (\ix ->
    --   let (r', c') = unlift $ unindex2 ix :: (Exp Int, Exp Int)
    --       r        = (r' `div` the sz) + (r' `mod` the sz)
    --       blockR   = r' `div` the sz
    --       blockC   = c' `div` the sz
    --       offset   = offsets ! lift (Z :. blockR :. blockC)
    --       c        = (c' * blockColCount) + (((c' `mod` blockColCount) + offset) `mod` the sz)
    --   in
    --   lift (Z :. r :. c)
    -- )
    -- arraylets

  -- undefined $
  -- A.filter (\p ->
  --   let (ix, v) = unlift p :: (Exp DIM2, Exp a)
  --   in
  --   undefined
  --   )
  --   (indexed arraylets)

  -- generate (lift (Z :. realRowCount)) $ \ix ->
  --   let r = unindex1 ix
  --   in
  --   sfoldl (\acc x ->
  --     f acc x
  --     )
  --     z
  --     ix
  --     arraylets

  -- permute
  --   (\p q ->
  --     let a = unlift p :: Exp a
  --         b = unlift q :: Exp a
  --     in
  --     f a b)
  --   (A.replicate (lift (Z :. realRowCount)) (unit z))
  --   (\ix ->
  --     let (r0, c0)         = unlift $ unindex2 ix :: (Exp Int, Exp Int)
  --         (blockR, blockC) = unlift $ unindex2 $ blockCoords !! c0 :: (Exp Int, Exp Int)
  --         offset           = offsets !! c0
  --     in index1 ((blockR*the sz) + r0)) $
  --   arraylets

  where
    (sz, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))
    (blockRowCount, blockColCount) = unlift $ unindex2 $ shape arrayletIndices :: (Exp Int, Exp Int)
    -- (blockRowCount, blockColCount) = unlift $ unindex2 $ shape arraylets :: (Exp Int, Exp Int) --(realRowCount `div` the sz, realColCount `div` the sz)

    -- (realRowCount, realColCount) = unlift $ unindex2 $ the realDim :: (Exp Int, Exp Int)
    realRowCount :: Exp Int
    realRowCount = blockRowCount*the sz

-- | A generalized self multiply-add.
matrixMultAdd :: forall a . (Elt a, a ~ Plain a) =>
  (Exp a -> Exp a -> Exp a) -> (Exp a -> Exp a -> Exp a) -> Exp a -> Exp a -> Acc (Matrixlet a) -> Acc (Scalar a)
matrixMultAdd mult add multZ addZ mLet =
  fold add addZ $
  fold mult multZ $

  generate (the realDim) 
    (\ix ->
      let (realR, realC)   = unlift $ unindex2 ix :: (Exp Int, Exp Int)
          (blockR, blockC) = (realR `div` sz, realC `div` sz)
          arrayletIx       = --(blockR*blockColCount) + blockC
            arrayletIndices ! lift (Z :. blockR :. blockC) :: Exp Int
          offset           = realR `mod` sz
      in
      cond (arrayletIx < 0)
           multZ
           (arraylets ! lift (Z :. arrayletIx :. offset))
    )

  -- backpermute
  --   (the realDim)
  --   (\ix ->
  --     let (realR, realC)   = unlift $ unindex2 ix :: (Exp Int, Exp Int)
  --         (blockR, blockC) = (realR `div` sz, realC `div` sz)
  --         arrayletIx       = (blockR*blockColCount) + blockC - unindex1 (zeroArrayletCounts ! lift (Z :. blockR :. blockC))
  --         offset           = realR `mod` sz
  --     in
  --     lift (Z :. arrayletIx :. offset)
  --   )

    -- (\ix ->
    --   let (r0, c0)         = unlift $ unindex2 ix :: (Exp Int, Exp Int)
    --       (blockR, blockC) = unlift $ unindex2 $ blockCoords !! c0 :: (Exp Int, Exp Int)
    --       offset           = offsets !! c0
    --   in index1 ((blockR*the sz) + r0)) $
    -- arraylets
  where
    sz = the sz0
    (blockRowCount, blockColCount) = unlift $ unindex2 $ shape arrayletIndices :: (Exp Int, Exp Int)
    (sz0, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))

-- -- TODO: Better name
-- foldColsMatrixlet' :: forall a . (Elt a, a ~ Plain a) => (Exp a -> Bool) -> (Exp a -> Exp a -> Exp a) -> Exp a -> Acc (Matrixlet a) -> Exp Bool
-- foldColsMatrixlet' continue f z mLet =
--   while
--     (continue . A.snd)
--     (\t ->
--       let (x, y) = unlift t :: (Exp a, Exp a)
--       in
--       f x y
--     )


matrixMatrixlet :: forall a b . (Elt a, Elt b) =>
  Acc (Matrixlet a) -> ((Exp Int, Exp Int) -> Exp b) -> Acc (Matrixlet b)
matrixMatrixlet mLet f =
  lift
    (sz0
    ,realDim
    ,blockCoords
    ,offsets
    ,arrayletIndices
    ,imap
       (\ix _ ->
          let (i, arrIx) = unlift $ unindex2 ix :: (Exp Int, Exp Int)
              (r0, c0)   = unlift $ unindex2 (blockCoords !! arrIx) :: (Exp Int, Exp Int)
              r          = (r0*sz) + i
              c          = (c0*sz) + ((i + (offsets !! arrIx)) `mod` sz)
          in
          f (r, c)
       )
       arraylets
     )
  where
    (sz0, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))

    sz = the sz0

    (realRowCount, realColCount) = unlift $ unindex2 $ the realDim :: (Exp Int, Exp Int)


mapMatrixlet :: forall a b . (Elt a, Elt b) =>
  (Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
mapMatrixlet f mLet =
  lift
    (sz0
    ,realDim
    ,blockCoords
    ,offsets
    ,arrayletIndices
    ,A.map f arraylets
    )
  where
    (sz0, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))

imapMatrixlet :: forall a b . (Elt a, Elt b) =>
  ((Exp Int, Exp Int) -> Exp a -> Exp b) -> Acc (Matrixlet a) -> Acc (Matrixlet b)
imapMatrixlet f mLet =
  lift
    (sz0
    ,realDim
    ,blockCoords
    ,offsets
    ,arrayletIndices
    -- ,imap (\ix v ->
    --         let (i, arrIx) = unlift $ unindex2 ix :: (Exp Int, Exp Int)
    --             (r0, c0)  = unlift $ unindex2 (blockCoords !! arrIx) :: (Exp Int, Exp Int)
    --             r         = (r0*sz) + i
    --             c         = (c0*sz) + ((i + (offsets !! arrIx)) `mod` sz)
    --         in
    --         f (r, c) v
    --       )
    --       arraylets
    -- )

    ,generate (shape arraylets)
      $ \ix ->
          let (i, arrIx) = unlift $ unindex2 ix :: (Exp Int, Exp Int)
              (r0, c0)  = unlift $ unindex2 (blockCoords !! arrIx) :: (Exp Int, Exp Int)
              r         = (r0*sz) + i
              c         = (c0*sz) + ((i + (offsets !! arrIx)) `mod` sz)
          in
          f (r, c) (arraylets ! ix)
    )
  where
    (sz0, realDim, blockCoords, offsets, arrayletIndices, arraylets) = unlift mLet
      :: (Acc (Scalar Int), Acc (Scalar DIM2), Acc (V DIM2), Acc (V Int), Acc (M Int), Acc (M a))

    sz = the sz0

    (realRowCount, realColCount) = unlift $ unindex2 $ the realDim :: (Exp Int, Exp Int)

ldpc' :: Acc (Matrixlet Double) -> Acc (Scalar Int) -> Acc (V Double) -> Acc (V Bool)
ldpc' mLet maxIterations orig_lam =
  map hard' $ loop 0 mLet orig_lam
  where
    loop :: Int -> Acc (Matrixlet Double) -> Acc (V Double) -> Acc (V Double)
    loop n ne lam =
      let (finalN, ne, r) = unlift (awhile loopCond loopBody liftedInit)
                            :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))
      in
      acond (the finalN >= the maxIterations)
            (lift orig_lam)
            r
      where
        liftedInit :: Acc (Scalar Int, Matrixlet Double, V Double)
        liftedInit = lift (unit (lift n), ne, lam)

    loopCond :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Bool)
    loopCond t = unit $ (the n < the maxIterations) && the (A.or ans) --(the done)
      where
        (n, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))

        -- done :: Acc (Scalar Bool)
        -- done =
        --   matrixMultAdd (/=) (||) (lift False) (lift False) $
        --   matrixMatrixlet mLet $ \ (r, c) -> hard' (lam !! c)

        ans :: Acc (V Bool)
        ans = foldColsMatrixlet (/=) (lift False) $ matrixMatrixlet mLet $ \ (r, c) -> hard' (lam !! c)


    loopBody :: Acc (Scalar Int, Matrixlet Double, V Double) -> Acc (Scalar Int, Matrixlet Double, V Double)
    loopBody t = lift (unit (the n+1), ne', lam')
      where
        (n, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (Matrixlet Double), Acc (V Double))

        ne_tanh'mat :: Acc (Matrixlet Double)
        ne_tanh'mat =
          imapMatrixlet
            (\ (m,n) v -> tanh (- ((lam !! n - v) / 2)))
            ne

        ne_tanhMulted :: Acc (V StableDiv)
        ne_tanhMulted = foldColsMatrixlet smult (lit 1) $ mapMatrixlet lit ne_tanh'mat

        ne' :: Acc (Matrixlet Double)
        ne' =
          imapMatrixlet
            (\ (m,n) v -> -2 * atanh'' ((ne_tanhMulted !! m) `sdiv` v))
            ne_tanh'mat

        lam' :: Acc (V Double)
        lam' = zipWith (+) orig_lam $ foldRowsMatrixlet (+) 0 ne'

