{-# LANGUAGE GADTs, DataKinds, KindSignatures, StandaloneDeriving, TypeFamilies, FlexibleInstances #-}

module Data.BitVector where

import Data.Bit

-- BitVector is a compact representation of a Shallow.
newtype BitVector = BitVector { unBitVector :: [Int] }
        deriving (Eq, Ord)

instance Show BitVector where
   show = concat . map show . toList

fromList :: [Bit] -> BitVector
fromList vs = BitVector [ i | (i,1) <- [1..] `zip` vs ]

toList :: BitVector -> [Bit]
toList (BitVector bv) = concat [ take (n-1) (repeat 0) ++ [1] | n <- zipWith (-) bv (0 : bv) ]

{-
--sshowBV :: BitVector -> String
--showBV (BitVector bv) = ""


zipWithBitVector :: (Bit -> Bit -> Bit) -> BitVector -> BitVector -> BitVector
zipWithBitVector f _ _ | f 0 0 /= 0 = error "zipWithBitVector assumes f 0 0 == 0"
zipWithBitVector f (BitVector xs) (BitVector ys) = BitVector [ z | (z,1) <- loop xs ys ]
  where
        loop []     []     = []
        loop (x:xs) []     = (x,f 1 0) : loop xs []
        loop []     (y:ys) = (y,f 0 1) : loop ys []
        loop (x:xs) (y:ys)
                | x == y   = (x,f 1 1) : loop xs ys
                | x < y    = (x,f 1 0) : loop xs (y:ys)
                | x > y    = (y,f 0 1) : loop (x:xs) ys
-}