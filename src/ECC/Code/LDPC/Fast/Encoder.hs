{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Fast.Encoder where

-- Fast implementations of LDPC encoder.

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Data.Monoid
import Debug.Trace

type M a = Matrix a
type V a = U.Vector a

encoder :: M Bool -> U.Vector Bool -> IO (U.Vector Bool)
encoder g v = return (v <> U.map toBool (U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) (fmap fromBool g)))))

--qc_encoder :: QuasiCyclic Word32 -> U.Vector Bool -> IO (U.Vector Bool)
--qc_encoder g v = return (v <> U.map toBool (U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) (fmap fromBool g)))))



