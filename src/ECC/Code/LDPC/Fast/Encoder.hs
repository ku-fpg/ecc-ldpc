{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Fast.Encoder where

-- Fast implementations of LDPC encoder.

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix as M
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Data.Monoid
import Debug.Trace
import qualified Data.Matrix.QuasiCyclic as Q
import Data.Word
import Data.Bits

type M a = Matrix a
type V a = U.Vector a

encoder :: Q.QuasiCyclic Integer -> Rate -> U.Vector Bool -> U.Vector Bool
encoder g@(Q.QuasiCyclic sz m) = 
    case sz of
      32 -> let m32 :: Matrix Word32
                m32 = fmap fromIntegral m 
            in \ r v -> 
                     let v' :: U.Vector Word32
                         v' = U.generate (nrows m) $ \ row -> 
                                 foldr (.|.) zeroBits [ if (v U.! (row * sz + i))
                                                        then bit i
                                                        else zeroBits
                                                      | i <- take sz [0..]
                                                      ]
                         res :: U.Vector Word32
                         res = U.generate (ncols m) $ \ col ->
                                     foldr xor zeroBits [ mulQWord32 (v' U.! row) (m32 M.! (row+1,col+1))
                                                        | row <- take (nrows m) [0..]
                                                        ]
                                 
                     in U.generate (sz * ncols m) $ \ col -> (res U.! (col `div` sz)) `testBit` (col `mod` sz)
      _  -> error $ "unsupported size for fast encoder : " ++ show sz


mulQWord32 :: Word32 -> Word32 -> Word32
mulQWord32 w1 w2 = go w2 0
  where
    go !w 32 = 0
    go !w !n | w1 `testBit` n = w `xor` go (w `rotateL` 1) (n+1)
             | otherwise      =         go (w `rotateL` 1) (n+1)

-- A reference implementation, to testing arbitrary finite bit width
refMulQWord32 :: Bits w => w -> w -> w
refMulQWord32 w1 w2 = foldr (.|.) zeroBits [ bit i | (i,True) <- take w [0..] `zip` U.toList r]
  where
    Just w = bitSizeMaybe w1
    g' = Q.toBitMatrix (Q.QuasiCyclic w $ M.fromLists [[w2]])
    v = U.generate w $ \ i -> w1 `testBit` i
    r = U.map toBool $ U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) ((fmap fromBool g'))))

