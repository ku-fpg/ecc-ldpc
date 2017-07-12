{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE RecordWildCards #-}

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
import Foreign.CUDA.Runtime.Marshal
import Foreign.CUDA.Types
-- import qualified Foreign.CUDA.Driver as CUDA
import Foreign.CUDA.Driver.Context.Base
import Foreign.CUDA.Driver.Exec
import qualified Foreign.CUDA.Driver.Device as CUDA
import Foreign.CUDA.Driver.Module

import GHC.Int

import Control.DeepSeq
import Data.Foldable (fold)
import Control.Monad (liftM2)

type IntT = Int32

data CudaAllocations =
  CudaAllocations
  { ldpcFun      :: Fun
  , cm           :: Module
  , ctx          :: Context
  -- , result_dev   :: DevicePtr IntT
  -- , offsets      :: DevicePtr IntT
  -- , mLet         :: DevicePtr Double
  -- , mLetRowCount :: IntT
  -- , mLetColCount :: IntT
  }

code :: Code
code = mkLDPC_CodeIO "cuda-arraylet-cm" E.encoder decoder initialize finalize

decoder ::
  CudaAllocations -> Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool))
decoder CudaAllocations{..} arr@(Q.QuasiCyclic sz _) rate maxIterations orig_lam  = do
  (mLet, offsets, rowCount, colCount, runningRowPop_dev, runningColPop_dev) <- init'd

  cm      <- loadFile "cudabits/cached_mult.ptx"
  ldpcFun <- getFun cm "ldpc"

  (orig_lam_dev, orig_lam_len) <- newListArrayLen $ U.toList $ orig_lam
  result_dev   <- mallocArray orig_lam_len :: IO (DevicePtr IntT)
  -- result_dev <- newListArray $ U.toList $ orig_lam

  launchKernel ldpcFun (orig_lam_len,1,1) (1,1,1) (orig_lam_len*8 + 16) Nothing [VArg mLet, IArg (fromIntegral sz), VArg runningRowPop_dev, VArg runningColPop_dev, VArg offsets, IArg rowCount, IArg colCount, IArg (fromIntegral maxIterations), VArg orig_lam_dev, IArg (fromIntegral orig_lam_len), VArg result_dev]

  sync

  result <- peekListArray orig_lam_len result_dev

  free result_dev
  free orig_lam_dev
  free mLet
  free offsets
  free runningRowPop_dev
  free runningColPop_dev

  -- unload cm

  return $! Just $! U.map toBool $! U.fromList result
  where
    init'd = initMatrixlet arr
    toBool 0 = False
    toBool 1 = True
    toBool n = error $ "toBool: invalid arg: " ++ show n

initialize :: IO CudaAllocations
initialize = do
  -- CUDA.initialise []
  -- dev0 <- device
  -- dev0props <- CUDA.props dev0
  -- print dev0props

  -- cm      <- loadFile "cudabits/cached_mult.ptx"
  -- ldpcFun <- getFun cm "ldpc"

  -- ctx     <- create dev0 []

  return (CudaAllocations {..})

finalize :: CudaAllocations -> IO ()
finalize CudaAllocations {..} = do
    return ()
  -- unload cm
  -- destroy ctx

initMatrixlet :: Q.QuasiCyclic Integer -> IO (DevicePtr Double, DevicePtr IntT, IntT, IntT, DevicePtr IntT, DevicePtr IntT)
initMatrixlet (Q.QuasiCyclic sz qm) = do
  mLetPtr    <- mallocArray (mLetRowCount * mLetColCount)
  offsetsPtr <- newListArray offsets

  memset mLetPtr (fromIntegral (mLetRowCount * mLetColCount)) 0

  let runningRowPop = map (fromIntegral :: Integer -> IntT) $ scanl (\x y -> x + pop y) 0 $ M.toLists qm
      runningColPop = map (fromIntegral :: Integer -> IntT) $ scanl (\x y -> x + pop y) 0 $ M.toLists $ M.transpose qm

  print (runningRowPop, runningColPop)

  runningRowPop_dev <- newListArray runningRowPop
  runningColPop_dev <- newListArray runningColPop

  return (mLetPtr, offsetsPtr, fromIntegral mLetRowCount, fromIntegral mLetColCount, runningRowPop_dev, runningColPop_dev)

  where
    mLetRowCount = sz
    mLetColCount = (M.nrows qm * M.ncols qm) - fromIntegral zeroArrayletCount

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


    zeroArrayletCount :: IntT
    zeroArrayletCount =
      getSum $
      foldMap (\n ->
        case n of
          0 -> 1
          _ -> 0) $
      qm

    offsets :: [IntT]
    offsets =
      foldMap (\n ->
        case n of
          0 -> []
          _ -> [g n]) $
      qm

