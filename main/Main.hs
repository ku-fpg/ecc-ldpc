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

codes :: Code
codes = Unencoded.code <> OrigReference.code <> Sparse.code <> Min.code <> SparseMin.code <> Model.code <> ElimTanh.code

-- usage: ./Main 0 2 4 6 8 0  bpsk
-- or, to run the LDPC reference implementation, at a single EBNO = 2.2
--        ./Main 2.2 ldpc/reference/jpl.1K/20
main :: IO ()
main = eccMain codes eccPrinter
