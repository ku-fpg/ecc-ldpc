{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RecordWildCards     #-}
{-# LANGUAGE BangPatterns        #-}

module ECC.Code.LDPC.GPU.CUDA.TwoArrays where

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
import           Foreign.CUDA.Runtime.Marshal as RM hiding (AllocFlag (..))
import Foreign.CUDA.Driver.Marshal (registerArray, AllocFlag (..))
import Foreign.CUDA.Types
-- import qualified Foreign.CUDA.Driver as CUDA
import Foreign.CUDA.Driver.Context.Base
import Foreign.CUDA.Driver.Exec
import qualified Foreign.CUDA.Driver.Device as CUDA
import qualified Foreign.CUDA.Driver.Stream as Stream
import Foreign.CUDA.Driver.Module

import Foreign.Storable (sizeOf)
import Foreign.ForeignPtr
import qualified Foreign.Marshal as F

import Data.IORef

import GHC.Int

import Control.DeepSeq
import Data.Foldable (fold)
import Control.Monad --(liftM2)

import GHC.Conc
import GHC.Float

type IntT    = Int32
type FloatTy = Float --Double

float_t_width :: Int
float_t_width = sizeOf (undefined :: FloatTy)

data CudaAllocations =
  CudaAllocations
  { cm           :: Module
  }

code :: Code
code = mkLDPC_CodeIO "two-arrays" 1 E.encoder decoder initialize finalize

pokeListArrayAsync :: S.Storable a => S.Vector a -> DevicePtr a -> Maybe Stream -> IO ()
pokeListArrayAsync !xs !dptr !stream = do
  let len       = S.length xs
      (fptr, _) = S.unsafeToForeignPtr0 xs
  withForeignPtr fptr (\p -> pokeArrayAsync len (HostPtr p) dptr stream)

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
  (mLet0, offsets, rowCount, colCount) <- init'd

  let rowsPerBlock
        | rowCount <= maxBlockSize = rowCount
        | otherwise                = rowCount `div` 2

      colsPerBlock
        | colCount <= maxBlockSize = colCount
        | otherwise                = colCount `div` 4

      rowBlockSize
        | rowCount <= maxBlockSize = rowCount `div` 2
        | otherwise                = rowCount `div` 4

      colBlockSize = colCount `div` 22

      productColsPerBlock = colBlockSize `div` 2

  makeNonzeroMatFun   <- getFun cm "makeNonzeroMat"
  tanhTransformFun    <- getFun cm "tanhTransform"
  setToOneFun         <- getFun cm "setToOne"
  insertOnesFun       <- getFun cm "insertOnes"
  selfProductFun      <- getFun cm "selfProduct"
  atanhTransformFun   <- getFun cm "atanhTransform"
  updateLamFun        <- getFun cm "updateLamT"
  parityRowResultsFun <- getFun cm "parityRowResults"
  checkParityFun      <- getFun cm "checkParity"

  memset mLet0 (fromIntegral (rowCount * colCount * fromIntegral float_t_width)) 0

  mLetRef <- newIORef mLet0

  newMLet0   <- mallocArray (fromIntegral $ rowCount * colCount) :: IO (DevicePtr FloatTy)
  newMLetRef <- newIORef newMLet0
  tempRef    <- newIORef =<< (mallocArray 1 :: IO (DevicePtr FloatTy))

  -- rowResults <- mallocArray (fromIntegral rowCount) :: IO (DevicePtr Int32)

  partials   <- mallocArray (fromIntegral rowCount) :: IO (DevicePtr FloatTy)
  print (rowCount, colCount)
  print (colCount*rowCount)

  mLetT   <- mallocArray (fromIntegral $ rowCount * colCount) :: IO (DevicePtr FloatTy)

  let orig_lam_len = fromIntegral colCount*sz :: Int

  orig_lam_dev <- mallocArray orig_lam_len :: IO (DevicePtr FloatTy)
  lam_dev <- mallocArray orig_lam_len :: IO (DevicePtr FloatTy)
  pop_dev <- newListArray [0] :: IO (DevicePtr Int32)
  done_dev <- newListArray [0] :: IO (DevicePtr Int32)

  stream1 <- Stream.create []
  stream2 <- Stream.create []
  mletStream <- Stream.create []

  -- putStr "atanh compute occupancy: "
  -- -- (fromIntegral (colCount `div` colBlockSize), fromIntegral (rowCount `div` rowBlockSize), 1)
  -- -- (fromIntegral colBlockSize,fromIntegral rowBlockSize,1)
  -- print (colBlockSize * rowBlockSize)

  return $ \rate maxIterations orig_lam -> do
    let orig_lam_stor = S.map double2Float $ U.convert orig_lam :: S.Vector FloatTy
        orig_lam_list = U.toList orig_lam

    -- pokeListArray orig_lam_list orig_lam_dev
    -- pokeListArray orig_lam_list lam_dev

    pokeListArrayAsync orig_lam_stor orig_lam_dev (Just stream1)
    pokeListArrayAsync orig_lam_stor lam_dev      (Just stream2)
    pokeListArray [0] pop_dev

    lamResultRef <- newIORef lam_dev

    mLet' <- readIORef mLetRef
    memset mLet' (fromIntegral $ rowCount * colCount * fromIntegral float_t_width) 0

    let go !iters
          | iters >= maxIterations = writeIORef lamResultRef orig_lam_dev
          | otherwise              = do
              mLet <- readIORef mLetRef
              newMLet <- readIORef newMLetRef

              -- Check parity every 5th iteration
              parity <- if True --iters < 3 || iters `rem` 5 == 0 || iters == maxIterations - 1
                then do
                  launchKernel parityRowResultsFun
                               (1, fromIntegral rowCount, 1)
                               (fromIntegral colCount, 1, 1)
                               0
                               Nothing
                               [VArg done_dev
                               ,VArg lam_dev
                               ,IArg rowCount
                               ,IArg colCount
                               ,IArg (fromIntegral sz)
                               ,VArg offsets
                               ]

                  [done] <- peekListArray 1 done_dev
                  return (done > 0)
                else return True

              writeIORef lamResultRef lam_dev

              when parity $ do
                -- Update matrix

                let nr = 4
                launchKernel selfProductFun
                             (fromIntegral rowCount `div` nr, 1, 1)
                             (nr, fromIntegral colCount, 1)
                             (fromIntegral colCount * nr * float_t_width)
                             Nothing
                             [VArg lam_dev
                             ,VArg mLet
                             ,VArg newMLet
                             ,VArg mLetT
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                copyArray orig_lam_len orig_lam_dev lam_dev

                swapRefs tempRef mLetRef newMLetRef

                -- Update guess

                launchKernel updateLamFun
                             -- (fromIntegral sz, 1, 1)
                             -- (fromIntegral colCount, 1, 1)
                             (fromIntegral (colCount `div` colBlockSize), fromIntegral (rowCount `div` rowBlockSize), 1)
                             (fromIntegral colBlockSize, fromIntegral rowBlockSize, 1)
                             -- (fromIntegral sz, 1, 1)
                             -- (fromIntegral colCount, 1, 1)
                             0
                             Nothing
                             [VArg lam_dev
                             ,VArg mLetT
                             ,IArg rowCount
                             ,IArg colCount
                             ,IArg (fromIntegral sz)
                             ,VArg offsets
                             ]

                go (iters+1)

    Stream.block stream1
    Stream.block stream2
    go 0

    lamPtr <- readIORef lamResultRef
    result <- peekListArray orig_lam_len lamPtr

    let r = Just $! U.map hard $! S.convert $! U.fromList result

    -- free orig_lam_dev
    -- free lam_dev
    -- free pop_dev

    -- enqueueFree orig_lam_dev
    -- enqueueFree lam_dev
    -- enqueueFree pop_dev

    return $! r
  where
    init'd = initMatrixlet arr

initialize :: IO CudaAllocations
initialize = do
  dummy <- RM.mallocArray 1 :: IO (DevicePtr Int) -- Sets up CUDA context
  cm    <- loadFile "cudabits/two_arrays.ptx"
  RM.free dummy

  return (CudaAllocations {..})

finalize :: CudaAllocations -> IO ()
finalize CudaAllocations {..} = do
    return ()

initMatrixlet :: Q.QuasiCyclic Integer -> IO (DevicePtr FloatTy, DevicePtr IntT, IntT, IntT)
initMatrixlet (Q.QuasiCyclic sz qm) = do
  mLetPtr    <- mallocArray (mLetRowCount * mLetColCount)
  offsetsPtr <- newListArray offsets

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
      foldMap (\n ->
        case n of
          0 -> [-1]
          _ -> [g n]) $
      qm

