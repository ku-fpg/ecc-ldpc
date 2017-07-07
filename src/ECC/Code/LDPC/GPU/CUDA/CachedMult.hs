{-# LANGUAGE ScopedTypeVariables #-}

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
import Foreign.CUDA.Driver.Exec -- (Fun(..), launchKernel)
import Foreign.CUDA.Driver.Module

import GHC.Int

code :: Code
code = mkLDPC_CodeIO "cuda-arraylet-cm" E.encoder decoder

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool))
decoder arr@(Q.QuasiCyclic sz _) = \rate maxIterations orig_lam -> do
  (mLet, offsets, rowCount, colCount) <- init'd

  cm   <- loadFile "cudabits/cached_mult.ptx"
  ldpc <- getFun cm "ldpc"

  result_dev   <- mallocArray (U.length orig_lam) :: IO (DevicePtr Int)
  orig_lam_dev <- newListArray $ U.toList $ orig_lam

  launchKernel ldpc (1,1,1) (1,1,1) 0 Nothing [VArg mLet, VArg offsets, IArg rowCount, IArg colCount, IArg (fromIntegral maxIterations), VArg orig_lam_dev, VArg result_dev]

  result <- peekListArray (U.length orig_lam) result_dev

  free mLet
  free offsets
  free result_dev
  free orig_lam_dev

  return $! Just $! U.map toBool $! U.fromList result
  where
    init'd = initMatrixlet arr
    toBool 0 = False
    toBool 1 = True
    toBool n = error $ "toBool: invalid arg: " ++ show n

initMatrixlet :: Q.QuasiCyclic Integer -> IO (DevicePtr Double, DevicePtr Int, Int32, Int32)
initMatrixlet (Q.QuasiCyclic sz qm) = do
  mLetPtr    <- mallocArray (mLetRowCount * mLetColCount)
  offsetsPtr <- newListArray offsets

  memset mLetPtr (fromIntegral (mLetRowCount * mLetColCount)) 0

  return (mLetPtr, offsetsPtr, fromIntegral mLetRowCount, fromIntegral mLetColCount)

  where
    mLetRowCount = sz
    mLetColCount = (M.nrows qm * M.ncols qm) - zeroArrayletCount

    -- The number must be a power of two, because there is only one bit set.
    g :: Integer -> Int
    g x | x `testBit` 0 = 0
        | x == 0        = error "got to zero; should never happen"
        | otherwise     = 1 + g (x `shiftR` 1)


    zeroArrayletCount :: Int
    zeroArrayletCount =
      getSum $
      foldMap (\n ->
        case n of
          0 -> 1
          _ -> 0) $
      qm

    offsets :: [Int]
    offsets =
      foldMap (\n ->
        case n of
          0 -> []
          _ -> [g n]) $
      qm

