{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}

module Data.Matrix.QuasiCyclic where

import Data.Bit
import Data.Bits
import Data.Matrix.Unboxed
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G
import Data.Char
import qualified Data.Matrix.Matlab as Matlab

import Data.Word

-- | Not abstract
--
-- You define a n-QuasiCyclic matrix by only including every n-th row.
--
data QuasiCyclic a = QuasiCyclic Int (Matrix a)

toBitMatrix :: (G.Vector U.Vector a, Bits a) => QuasiCyclic a -> Matrix Bool
toBitMatrix (QuasiCyclic sz a) = fromLists
        [ [ (a ! (m,n)) `testBit` ((j - i) `mod` sz)
          | n <- [1..cols a], j <- [0..sz-1]
          ]
        | m <- [1..rows a], i <- [0..sz-1]
        ]

fromBitMatrix :: (G.Vector U.Vector a, Bits a) => Int -> Matrix Bool -> QuasiCyclic a
fromBitMatrix sz a | toBitMatrix q == a = q
                   | otherwise          = error "fromBitMatrix: not a cyclic matrix"
        -- We take every n'th line, generate a QuasiCyclic matrix,
        -- and then test the input against the expanded candidate result
   where
          q = QuasiCyclic sz $ fromLists
                [ [ foldr (.|.) zero
                    [ bit j
                    | j <- [0..sz-1]
                    , let v = a ! (m,n + j)
                    , v  -- if set
                    ]
                  | n <- [1,(sz+1)..cols a]
                  ]
                | m <- [1,(sz+1)..rows a]
                ]
          zero = bit 0 `xor` bit 0 -- Hack to get a zero, without


cycleSize :: Bits a => QuasiCyclic a -> Int
cycleSize (QuasiCyclic n a) = n

instance (G.Vector U.Vector a, Show a) => Show (QuasiCyclic a) where
        -- The first line is the cycle size; the rest is a Matlab style matrix
        show (QuasiCyclic sz a) = show sz ++ "\n" ++ show (Matlab.Matlab a)

instance (G.Vector U.Vector a, Read a) => Read (QuasiCyclic a) where
        readsPrec _ txt0 = [ (QuasiCyclic sz a,txt2)
                           | (sz,txt1)              <- reads (dropWhile isSpace txt0)
                           , (Matlab.Matlab a,txt2) <- reads (dropWhile isSpace txt1)
                           ]


example :: QuasiCyclic Word8
example = QuasiCyclic 8 $ fromLists [[ 0x01, 0xf8], [0x0, 0xd1]]

bitM :: Matrix Bool
bitM = fromLists [[False,True,False,True,True,False],[False,False,True,False,True,True],[True,False,False,True,False,True]]
