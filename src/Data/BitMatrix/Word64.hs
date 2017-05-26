module Data.BitMatrix.Word64 where

import Data.Bit
import Data.Bits
import Data.BitVector.Word64 as BV

import qualified Data.Vector as V

-- The matrix is indexed start at (1,1)
newtype BitMatrix = BitMatrix { unBitMatrix :: V.Vector BV.BitVector }
        deriving Show

fromLists :: [[Bool]] -> BitMatrix
fromLists = BitMatrix . V.fromList . map BV.fromList

toLists :: BitMatrix -> [[Bool]]
toLists = fmap BV.toList . V.toList . unBitMatrix

(!) :: BitMatrix -> (Int,Int) -> Bool
(!) (BitMatrix vs) (m,n) = (vs V.! (m - 1)) BV.! (succ (n - 1))

vecMatMul :: BitVector -> BitMatrix -> BitVector
vecMatMul bv bm = foldr1 (zipWithWord64 xor) [ row | (True,row) <- zip (BV.toList bv) (V.toList (unBitMatrix bm)) ]

-- intentually lazy in result
parityMatVecMul :: BitMatrix -> BitVector -> [Bool]
parityMatVecMul bm bv = [ odd $ BV.popCount $ zipWithWord64 (.&.) bv row | row <- V.toList (unBitMatrix bm) ]


m = fromLists [[True,False,False,True,False,True]
              ,[False,True,False,True,True,True]
              ,[False,False,True,True,True,False]]

v = BV.fromList [True,False,True,False,True,True]

