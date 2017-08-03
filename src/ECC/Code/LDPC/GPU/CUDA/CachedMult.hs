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
code = mkLDPC_CodeIO "cuda-arraylet-cm" E.encoder decoder initialize finalize


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

decoder ::
  CudaAllocations -> Q.QuasiCyclic Integer -> IO (Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool)))
decoder CudaAllocations{..} arr@(Q.QuasiCyclic sz _) = do
  (mLet0, offsets, rowCount, colCount) <- init'd

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

  -- nonzeroMat <- mallocArray (fromIntegral $ rowCount * colCount) :: IO (DevicePtr Bool)
  partials   <- mallocArray (fromIntegral rowCount) :: IO (DevicePtr FloatTy)
  print (rowCount, colCount)
  print (colCount*rowCount)

  mLetT   <- mallocArray (fromIntegral $ rowCount * colCount) :: IO (DevicePtr FloatTy)

  -- launchKernel makeNonzeroMatFun
  --              (fromIntegral colCount, 1, 1)
  --              (1, fromIntegral rowCount, 1)
  --              8
  --              Nothing
  --              [VArg nonzeroMat
  --              ,IArg rowCount
  --              ,IArg colCount
  --              ,IArg (fromIntegral sz)
  --              ,VArg offsets
  --              ]


  return $ \rate maxIterations orig_lam -> do
    let orig_lam_list = convertToFloatT $ U.toList orig_lam :: [FloatTy]
    (orig_lam_dev, orig_lam_len) <- newListArrayLen orig_lam_list

    lam_dev <- newListArray orig_lam_list

    pop_dev <- newListArray [0] :: IO (DevicePtr Int32)

    lamResultRef <- newIORef lam_dev

    mLet' <- readIORef mLetRef
    memset mLet' (fromIntegral $ rowCount * colCount * fromIntegral float_t_width) 0

    stream1 <- Stream.create []

    let go !iters
          | iters >= maxIterations = writeIORef lamResultRef orig_lam_dev
          | otherwise              = do
              mLet <- readIORef mLetRef
              newMLet <- readIORef newMLetRef

              -- Check parity
              launchKernel parityRowResultsFun
                           (1, fromIntegral rowCount, 1)
                           (fromIntegral colCount, 1, 1)
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
                           (1,1,1)
                           (1, fromIntegral rowCount, 1)
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
                launchKernel setToOneFun
                             (fromIntegral colCount, 1, 1)
                             (1, fromIntegral rowCount, 1)
                             0
                             (Just stream1)
                             [VArg newMLet
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]
                launchKernel tanhTransformFun
                             (fromIntegral colCount, 1, 1)
                             (1, fromIntegral rowCount, 1)
                             0
                             Nothing
                             [VArg mLet
                             ,VArg lam_dev
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]
                Stream.block stream1

                -- launchKernel insertOnesFun
                --              (fromIntegral colCount, 1, 1)
                --              (1, fromIntegral rowCount, 1)
                --              0
                --              Nothing
                --              [VArg mLet
                --              -- ,VArg nonzeroMat
                --              ,IArg rowCount
                --              ,IArg colCount
                --              ,IArg (fromIntegral sz)
                --              ,VArg offsets
                --              ]

                -- NOTE: Assumes column count is divisible by 11 and row
                -- count divisible by 2.

                launchKernel selfProductFun
                             (1, 2, fromIntegral colCount)
                             -- (fromIntegral (colCount `div` 4), fromIntegral (rowCount `div` 8), 1)
                             (fromIntegral (colCount `div` 22), fromIntegral (rowCount `div` 2), 1)
                             -- (1, fromIntegral rowCount, fromIntegral colCount)
                             -- (fromIntegral colCount, 1, 1)
                             (fromIntegral colCount * float_t_width)
                             Nothing
                             [VArg mLet
                             ,VArg newMLet
                             -- ,VArg nonzeroMat
                             -- ,VArg lam_dev
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                -- launchKernel selfProductRowsFun
                --              (fromIntegral rowCount, 1, 1)
                --              (fromIntegral colCount, 1, 1)
                --              ((fromIntegral colCount+1)*8)
                --              Nothing
                --              [VArg mLet
                --              ,VArg newMLet
                --              -- ,VArg nonzeroMat
                --              -- ,VArg lam_dev
                --              ,IArg rowCount
                --              ,IArg colCount
                --              ,IArg (fromIntegral sz)
                --              ,VArg offsets
                --              ]

                launchKernel atanhTransformFun
                             (11, 2, 1)
                             (fromIntegral (colCount`div`11),fromIntegral (rowCount `div` 2),1)
                             -- (fromIntegral colCount, 1, 1)
                             -- (1,fromIntegral rowCount,1)
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
                             (11, 2, 1)
                             (fromIntegral (colCount`div`11),fromIntegral (rowCount `div` 2),1)
                             -- (fromIntegral colCount, 1, 1)
                             -- (1, fromIntegral rowCount, 1)
                             --
                             -- (1, fromIntegral rowCountT, 1)
                             -- (fromIntegral colCountT,1,1)
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

initMatrixlet :: Q.QuasiCyclic Integer -> IO (DevicePtr FloatTy, DevicePtr IntT, IntT, IntT)
initMatrixlet (Q.QuasiCyclic sz qm) = do
  mLetPtr    <- mallocArray (mLetRowCount * mLetColCount)
  -- mLetPtr    <- mallocManagedArray [CuMemAttachHost] (mLetRowCount * mLetColCount)
  offsetsPtr <- newListArray offsets

  -- memset mLetPtr (fromIntegral (mLetRowCount * mLetColCount * 8)) 0

  return (mLetPtr, offsetsPtr, fromIntegral mLetRowCount, fromIntegral mLetColCount)

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


    offsets :: [IntT]
    offsets =
      -- traceShowId $
      foldMap (\n ->
        case n of
          0 -> [-1]
          -- 0 -> [0]
          _ -> [g n]) $
      qm

