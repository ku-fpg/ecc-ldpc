module Main where

import ECC.Tester
import ECC.Types
import Data.Monoid

import qualified ECC.Code.BPSK as BPSK
import qualified ECC.Code.LDPC.Reference as Reference

codes :: Code
codes = BPSK.code <> Reference.code

-- usage: ./Main 0 2 4 6 8 0  bpsk
main :: IO ()
main = eccMain codes eccPrinter
