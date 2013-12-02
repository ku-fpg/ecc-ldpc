{-# LANGUAGE GADTs, DataKinds, KindSignatures, StandaloneDeriving, TypeFamilies, FlexibleInstances #-}

module Data.BitVector.Sparse where

import Data.Bit
import qualified Data.List as L

-- BitVector is a compact representation of a Sparse bit vector.
-- It has not specific length, other than where the 1's are.
newtype BitVector = BitVector { unBitVector :: [Int] }
        deriving (Eq, Ord)

instance Show BitVector where
   show = concat . map show . toList

fromList :: [Bit] -> BitVector
fromList vs = BitVector [ i | (i,1) <- [1..] `zip` vs ]

toList :: BitVector -> [Bit]
toList (BitVector bv) = concat [ take (n-1) (repeat 0) ++ [1] | n <- L.zipWith (-) bv (0 : bv) ]

elems :: BitVector -> [Int]
elems = unBitVector

vector :: [Int] -> BitVector
vector = BitVector . L.sort


-- zipWith combines two bit vectors with a binary bit function
zipWith :: (Bit -> Bit -> Bit) -> BitVector -> BitVector -> BitVector
zipWith f _ _ | f 0 0 /= 0 = error "zipWith assumes f 0 0 == 0"
zipWith f (BitVector xs) (BitVector ys) = BitVector [ z | (z,1) <- loop xs ys ]
  where
        loop []     []     = []
        loop (x:xs) []     = (x,f 1 0) : loop xs []
        loop []     (y:ys) = (y,f 0 1) : loop ys []
        loop (x:xs) (y:ys)
                | x == y   = (x,f 1 1) : loop xs ys
                | x < y    = (x,f 1 0) : loop xs (y:ys)
                | x > y    = (y,f 0 1) : loop (x:xs) ys


(!) :: BitVector -> Int -> Bit
(!) (BitVector xs) n = if n `elem` xs then 1 else 0
