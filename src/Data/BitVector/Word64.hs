module Data.BitVector.Word64 where

import Data.Bits as B
import Data.Bit
import Data.List
import Data.Word

import qualified Data.Vector.Unboxed as U

-- Perhaps this should be a vector?
newtype BitVector = BitVector { unBitVector :: U.Vector Word64 }
        deriving Show

fromList :: [Bit] -> BitVector
fromList = BitVector
         . U.fromList
         . unfoldr (\ ws -> if null ws then Nothing else Just (parseBits ws))


toList :: BitVector -> [Bit]
toList (BitVector ws) = concat [ [ mkBit $ testBit w i | i <- [0..63] ] | w <- U.toList ws ]


parseBits :: (Num b, Bits b) => [Bit] -> (b,[Bit])
parseBits bs = (r,drop sz bs)
  where r = foldr (.|.) 0 [ setBit 0 i | (i,1) <- zip [0..sz-1] bs ]
        sz = bitSize r

(!) :: BitVector -> Int -> Bit
(!) (BitVector vs) n = mkBit $ testBit (vs U.! (n `div` 64)) (n `mod` 64)

-- zipWith combines two bit vectors with a binary bit function
zipWithWord64 :: (Word64 -> Word64 -> Word64) -> BitVector -> BitVector -> BitVector
zipWithWord64 f (BitVector bv1) (BitVector bv2) = BitVector (U.zipWith f bv1 bv2)

popCount :: BitVector -> Int
popCount = sum . map B.popCount . U.toList . unBitVector

