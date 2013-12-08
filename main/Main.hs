module Main where

import ECC.Tester
import ECC.Types
import Data.Monoid

import qualified ECC.Code.BPSK as BPSK
import qualified ECC.Code.LDPC.Reference as Reference
import qualified ECC.Code.LDPC.Model as Model

import qualified ECC.Code.LDPC.ElimTanh as ElimTanh

codes :: Code
codes = BPSK.code <> Reference.code <> Model.code <> ElimTanh.code

-- usage: ./Main 0 2 4 6 8 0  bpsk
-- or, to run the LDPC reference implementation, at a single EBNO = 2.2
--        ./Main 2.2 ldpc/reference/jpl.1K/20
main :: IO ()
main = eccMain codes eccPrinter
