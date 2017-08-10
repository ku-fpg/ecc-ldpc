{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, TypeApplications #-}
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
import Data.DoubleWord
import Data.Bits
import Data.Proxy

type M a = Matrix a
type V a = U.Vector a

encoder :: Q.QuasiCyclic Integer -> Rate -> U.Vector Bool -> U.Vector Bool
encoder g@(Q.QuasiCyclic sz m) =
    case sz of
      32  -> go (Proxy @Word32)
      64  -> go (Proxy @Word64)
      128 -> go (Proxy @Word128)
      256 -> go (Proxy @Word256)
      _   -> error $ "unsupported size for fast encoder : " ++ show sz
  where
    go :: forall w. (Num w, FiniteBits w) =>
      Proxy w -> Rate -> U.Vector Bool -> U.Vector Bool
    go Proxy =
      let m32 :: Matrix w
          m32 = fmap fromIntegral m
      in \ r v ->
               let v' :: V.Vector w
                   v' = V.generate (nrows m) $ \ row ->
                           foldr (.|.) zeroBits [ if (v U.! (row * sz + i))
                                                  then bit i
                                                  else zeroBits
                                                | i <- take sz [0..]
                                                ]
                   res :: V.Vector w
                   res = V.generate (ncols m) $ \ col ->
                               foldr xor zeroBits [ mulWord (v' V.! row) (m32 M.! (row+1,col+1))
                                                  | row <- take (nrows m) [0..]
                                                  ]

               in
               V.convert $ V.generate (sz * ncols m) $ \ col -> (res V.! (col `div` sz)) `testBit` (col `mod` sz)


mulWord :: (Num a, FiniteBits a) => a -> a -> a
mulWord w1 w2 = go w2 0
  where
    go !w !n | n == finiteBitSize w1 = 0
    go !w !n | w1 `testBit` n        = w `xor` go (w `rotateL` 1) (n+1)
             | otherwise             = go (w `rotateL` 1) (n+1)

-- A reference implementation, to testing arbitrary finite bit width
refMulQWord32 :: Bits w => w -> w -> w
refMulQWord32 w1 w2 = foldr (.|.) zeroBits [ bit i | (i,True) <- take w [0..] `zip` U.toList r]
  where
    Just w = bitSizeMaybe w1
    g' = Q.toBitMatrix (Q.QuasiCyclic w $ M.fromLists [[w2]])
    v = U.generate w $ \ i -> w1 `testBit` i
    r = U.map toBool $ U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) ((fmap fromBool g'))))

