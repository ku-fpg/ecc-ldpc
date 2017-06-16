{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, DeriveFunctor #-}
module ECC.Code.LDPC.Zero where

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Data.Bit
import Data.Bits
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix.QuasiCyclic as Q

import qualified ECC.Code.LDPC.Fast.Encoder as E

type M a = Matrix a
type V a = U.Vector a

code :: Code
code = mkLDPC_Code "ldpc-zero" E.encoder decoder

---------------------------------------------------------------------

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a = \ _rate _maxIterations orig_lam -> Just (U.map hard orig_lam)

