{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RecordWildCards #-}
{-# LANGUAGE BangPatterns #-}

module ECC.Code.LDPC.GPU.CUDA.CachedMult where

-- Uses the StableDiv data structure to cache multiplications.
-- Uses CUDA directly to run on the GPU

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import qualified Data.Matrix as M
import Data.Bit
import Data.Bits
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Storable as S
import qualified Data.Matrix.QuasiCyclic as Q
import Debug.Trace

import qualified ECC.Code.LDPC.Fast.Encoder as E

import Data.Monoid

-- import Foreign.CUDA hiding (launchKernel)
import           Foreign.CUDA.Runtime.Marshal as RM
import Foreign.CUDA.Types
-- import qualified Foreign.CUDA.Driver as CUDA
import Foreign.CUDA.Driver.Context.Base
import Foreign.CUDA.Driver.Exec
import qualified Foreign.CUDA.Driver.Device as CUDA
import qualified Foreign.CUDA.Driver.Stream as Stream
import Foreign.CUDA.Driver.Module

import Foreign.Storable (sizeOf)

import Data.IORef

import GHC.Int

import Control.DeepSeq
import Data.Foldable (fold)
import Control.Monad --(liftM2)

import Data.List (sortBy)
import Data.Ord

type IntT    = Int32
type FloatTy = Double

float_t_width :: Int
float_t_width = sizeOf (undefined :: FloatTy)

data CudaAllocations =
  CudaAllocations
  { cm           :: Module
  -- , result_dev   :: DevicePtr IntT
  -- , offsets      :: DevicePtr IntT
  -- , mLet         :: DevicePtr FloatTy
  -- , rowCount :: IntT
  -- , colCount :: IntT
  }

code :: Code
code = mkLDPC_CodeIO "cuda-arraylet-cm" 1 E.encoder decoder initialize finalize


                -- atomicWriteIORef tempRef mLet
                -- atomicWriteIORef mLetRef newMLet
                -- tempPtr <- readIORef tempRef
                -- atomicWriteIORef newMLetRef tempPtr

swapRefs :: IORef a -> IORef a -> IORef a -> IO ()
swapRefs tempRef xRef yRef = do
  x <- readIORef xRef
  y <- readIORef yRef

  atomicWriteIORef tempRef x
  atomicWriteIORef xRef y

  temp <- readIORef tempRef

  atomicWriteIORef yRef temp


convertToFloatT :: [Double] -> [FloatTy]
convertToFloatT = map realToFrac

maxBlockSize :: IntT
maxBlockSize = 1024

decoder ::
  CudaAllocations -> Q.QuasiCyclic Integer -> IO (Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool)))
decoder CudaAllocations{..} arr@(Q.QuasiCyclic sz _) = do
  (mLet0, offsets, nonzeroIndicesR, nonzeroIndicesC, nonzeroIndicesLen, rowCount, colCount) <- init'd

  let rowsPerBlock
        | rowCount <= maxBlockSize = rowCount
        | otherwise                = rowCount `div` 2

      colsPerBlock
        | colCount <= maxBlockSize = colCount
        | otherwise                = colCount `div` 4

      rowBlockSize
        | rowCount <= maxBlockSize = rowCount `div` 2
        | otherwise                = rowCount `div` 8

      colBlockSize = colCount `div` 11

      productColsPerBlock = colBlockSize `div` 2

  makeNonzeroMatFun <- getFun cm "makeNonzeroMat"
  tanhTransformFun  <- getFun cm "tanhTransform"
  setToOneFun       <- getFun cm "setToOne"
  insertOnesFun     <- getFun cm "insertOnes"
  selfProductFun    <- getFun cm "selfProduct"
  selfProductRowsFun <- getFun cm "selfProductRows"
  atanhTransformFun <- getFun cm "atanhTransform"
  updateLamFun      <- getFun cm "updateLam"
  parityRowResultsFun <- getFun cm "parityRowResults"
  checkParityFun    <- getFun cm "checkParity"

  memset mLet0 (fromIntegral (rowCount * colCount * fromIntegral float_t_width)) 0

  mLetRef <- newIORef mLet0

  newMLet0   <- mallocArray (fromIntegral $ rowCount * colCount) :: IO (DevicePtr FloatTy)
  newMLetRef <- newIORef newMLet0
  tempRef    <- newIORef =<< (mallocArray 1 :: IO (DevicePtr FloatTy))

  rowResults <- mallocArray (fromIntegral rowCount) :: IO (DevicePtr Bool)

  print (rowCount, colCount)
  print (colCount*rowCount)

  return $ \rate maxIterations orig_lam -> do
    let orig_lam_list = convertToFloatT $ U.toList orig_lam :: [FloatTy]
    (orig_lam_dev, orig_lam_len) <- newListArrayLen orig_lam_list

    lam_dev <- newListArray orig_lam_list

    pop_dev <- newListArray [0] :: IO (DevicePtr Int32)

    lamResultRef <- newIORef lam_dev

    mLet' <- readIORef mLetRef
    memset mLet' (fromIntegral $ rowCount * colCount * fromIntegral float_t_width) 0

    let go !iters
          | iters >= maxIterations = writeIORef lamResultRef orig_lam_dev
          | otherwise              = do
              mLet <- readIORef mLetRef
              newMLet <- readIORef newMLetRef

              -- Check parity
              launchKernel parityRowResultsFun
                           (fromIntegral (colCount `div` colsPerBlock), fromIntegral rowCount, 1)
                           (fromIntegral colsPerBlock, 1, 1)
                           -- (1, fromIntegral rowCount, 1)
                           -- (fromIntegral colCount, 1, 1)
                           float_t_width
                           Nothing
                           [VArg rowResults
                           ,VArg lam_dev
                           ,IArg rowCount
                           ,IArg colCount
                           ,IArg (fromIntegral sz)
                           ,VArg offsets
                           ]

              launchKernel checkParityFun
                           (1, fromIntegral (rowCount `div` rowsPerBlock),1)
                           (1, fromIntegral rowsPerBlock, 1)
                           float_t_width
                           Nothing
                           [VArg pop_dev
                           ,VArg rowResults
                           ]

              [pop] <- peekListArray 1 pop_dev
              let parity = pop > 0

              writeIORef lamResultRef lam_dev

              when parity $ do
                -- Update matrix
                launchKernel tanhTransformFun
                             (fromIntegral colCount, fromIntegral (rowCount `div` rowsPerBlock), 1)
                             (1, fromIntegral rowsPerBlock, 1)
                             0
                             Nothing
                             [VArg mLet
                             ,VArg lam_dev
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                -- let rowsPerBlock = 4
                let indicesPerBlock = 1
                launchKernel selfProductFun
                             -- (1, fromIntegral rowCount `div` rowsPerBlock, fromIntegral colCount)
                             (fromIntegral nonzeroIndicesLen `div` indicesPerBlock, 1, 1)
                             (fromIntegral colCount, indicesPerBlock, 1)
                             (fromIntegral colCount * indicesPerBlock * float_t_width)
                             Nothing
                             [VArg mLet
                             ,VArg newMLet
                             ,VArg nonzeroIndicesR
                             ,VArg nonzeroIndicesC
                             ,IArg nonzeroIndicesLen
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                launchKernel atanhTransformFun
                             (fromIntegral (colCount `div` colBlockSize), fromIntegral (rowCount `div` rowBlockSize), 1)
                             (fromIntegral colBlockSize,fromIntegral rowBlockSize,1)
                             0
                             Nothing
                             [VArg newMLet
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                swapRefs tempRef mLetRef newMLetRef

                -- Update guess
                copyArray orig_lam_len orig_lam_dev lam_dev

                launchKernel updateLamFun
                             (fromIntegral (colCount `div` colBlockSize), fromIntegral (rowCount `div` rowBlockSize), 1)
                             (fromIntegral colBlockSize,fromIntegral rowBlockSize,1)
                             0
                             Nothing
                             [VArg lam_dev
                             ,VArg newMLet
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                go (iters+1)

    go 0

    lamPtr <- readIORef lamResultRef
    result <- peekListArray orig_lam_len lamPtr

    let r = Just $! U.map hard $! S.convert $! U.fromList result

    free orig_lam_dev

    -- free mLet
    -- free newMLet
    -- free offsets
    -- free zeroArr
    -- free pop_dev

    -- unload cm

    return $! r
  where
    init'd = initMatrixlet arr

initialize :: IO CudaAllocations
initialize = do
  dummy <- RM.mallocArray 1 :: IO (DevicePtr Int) -- Sets up CUDA context
  cm    <- loadFile "cudabits/cached_mult.ptx"
  RM.free dummy

  return (CudaAllocations {..})

finalize :: CudaAllocations -> IO ()
finalize CudaAllocations {..} = do
    return ()
  -- unload cm
  -- destroy ctx

initMatrixlet :: Q.QuasiCyclic Integer -> IO (DevicePtr FloatTy, DevicePtr IntT, DevicePtr IntT, DevicePtr IntT, IntT, IntT, IntT)
initMatrixlet (Q.QuasiCyclic sz qm) = do
  mLetPtr    <- mallocArray (mLetRowCount * mLetColCount)
  -- mLetPtr    <- mallocManagedArray [CuMemAttachHost] (mLetRowCount * mLetColCount)
  offsetsPtr <- newListArray offsets
  (nonzeroIndicesRPtr, nonzeroIndicesLen) <- newListArrayLen nonzeroIndicesR
  nonzeroIndicesCPtr <- newListArray nonzeroIndicesC

  putStrLn $ "nonzeroIndicesLen = " ++ show nonzeroIndicesLen
  putStrLn $ "total index count = " ++ show (length offsets*sz)

  -- memset mLetPtr (fromIntegral (mLetRowCount * mLetColCount * 8)) 0

  return (mLetPtr, offsetsPtr, nonzeroIndicesRPtr, nonzeroIndicesCPtr, fromIntegral nonzeroIndicesLen, fromIntegral mLetRowCount, fromIntegral mLetColCount)

  where
    mLetRowCount = M.nrows qm*sz
    mLetColCount = M.ncols qm

    pop :: [Integer] -> Integer
    pop =
      getSum .
      foldMap (\n ->
        if n == 0
        then 0
        else 1)

    -- The number must be a power of two, because there is only one bit set.
    g :: Integer -> IntT
    g x | x `testBit` 0 = 0
        | x == 0        = error "got to zero; should never happen"
        | otherwise     = 1 + g (x `shiftR` 1)

    nonzeroIndicesR :: [IntT]
    nonzeroIndicesR = map fst nonzeroIndices

    nonzeroIndicesC :: [IntT]
    nonzeroIndicesC = map snd nonzeroIndices

    nonzeroIndices :: [(IntT, IntT)]
    nonzeroIndices =
      sortBy (\(r1, c1) (r2, c2) -> compare (c1, r1) (c2, r2)) $
      foldMap (\ (r, c) ->
        let v   = qm M.! (fromIntegral $ (r `div` sz)+1, fromIntegral $ c+1)
            off = g v
        in
        if v /= 0
        then [(fromIntegral r, fromIntegral c)]
        -- then [(fromIntegral $ r+i, fromIntegral $ (c + ((i + off) `mod` fromIntegral sz) `mod` fromIntegral sz)) | i <- [0..fromIntegral sz-1]]
        else [])
      [ (r, c) | r <- [0..fromIntegral mLetRowCount-1], c <- [0..fromIntegral mLetColCount-1]]

    offsets :: [IntT]
    offsets =
      -- traceShowId $
      foldMap (\n ->
        case n of
          0 -> [-1]
          -- 0 -> [0]
          _ -> [g n]) $
      qm

