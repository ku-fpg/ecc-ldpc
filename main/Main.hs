{-# LANGUAGE OverloadedStrings #-}

module Main where

import ECC.Tester
import ECC.Types
import Data.Monoid

import qualified ECC.Code.Unencoded as Unencoded
import qualified ECC.Code.LDPC.Reference.Orig as OrigReference
import qualified ECC.Code.LDPC.Reference.Sparse as Sparse
import qualified ECC.Code.LDPC.Reference.Min as Min
import qualified ECC.Code.LDPC.Reference.SparseMin as SparseMin
import qualified ECC.Code.LDPC.Model as Model

import qualified ECC.Code.LDPC.ElimTanh as ElimTanh
import qualified ECC.Code.LDPC.Fast.Arraylet as Arraylet
import qualified ECC.Code.LDPC.Fast.ArrayletMin as ArrayletMin
import qualified ECC.Code.LDPC.Fast.CachedMult as CachedMult
import qualified ECC.Code.LDPC.GPU.CachedMult as GPUCachedMult
import qualified ECC.Code.LDPC.GPU.Reference as GPUReference
import qualified ECC.Code.LDPC.Zero as Zero
import qualified ECC.Code.LDPC.GPU.CUDA.CachedMult as CUDACachedMult
import qualified ECC.Code.LDPC.GPU.CUDA.Arraylet1 as CUDAArraylet1
import qualified ECC.Code.LDPC.GPU.CUDA.Arraylet2 as CUDAArraylet2
import qualified ECC.Code.LDPC.GPU.CUDA.TwoArrays as TwoArrays

import Data.Array.Accelerate.LLVM.PTX

import Data.Array.Accelerate.Debug
import System.Metrics
import System.Remote.Monitoring

codes :: Code
codes = Unencoded.code <> OrigReference.code <> Sparse.code <> Min.code <> SparseMin.code <> Model.code <> ElimTanh.code
            <> Arraylet.code <> ArrayletMin.code <> CachedMult.code <> GPUCachedMult.code <> GPUReference.code <> Zero.code <> CUDACachedMult.code <> CUDAArraylet1.code <> CUDAArraylet2.code <> TwoArrays.code

-- usage: ./Main 0 2 4 6 8 0  bpsk
-- or, to run the LDPC reference implementation, at a single EBNO = 2.2
--        ./Main 2.2 ldpc/reference/jpl.1K/20
main :: IO ()
main = do
  -- store <- initAccMetrics
  -- registerGcMetrics store
  -- server <- forkServerWith store "localhost" 8000

  -- registerPinnedAllocator
  eccMain codes eccPrinter
