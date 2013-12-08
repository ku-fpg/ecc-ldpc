module Data.BitMatrix.Word64 where

import Data.Bit
import Data.Bits
import Data.BitVector.Word64 as BV

import qualified Data.Vector as V

-- The matrix is indexed start at (1,1)
newtype BitMatrix = BitMatrix { unBitMatrix :: V.Vector BV.BitVector }
        deriving Show

fromLists :: [[Bit]] -> BitMatrix
fromLists = BitMatrix . V.fromList . map BV.fromList

toLists :: BitMatrix -> [[Bit]]
toLists = fmap BV.toList . V.toList . unBitMatrix

(!) :: BitMatrix -> (Int,Int) -> Bit
(!) (BitMatrix vs) (m,n) = (vs V.! (m - 1)) BV.! (succ (n - 1))

vecMatMul :: BitVector -> BitMatrix -> BitVector
vecMatMul bv bm = foldr1 (zipWithWord64 xor) [ row | (1,row) <- zip (BV.toList bv) (V.toList (unBitMatrix bm)) ]

m = fromLists [[1,0,0,1,0,1]
              ,[0,1,0,1,1,1]
              ,[0,0,1,1,1,0]]

v = BV.fromList [1,0,1,0,1,1]

