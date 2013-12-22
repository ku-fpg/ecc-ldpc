module Data.Matrix.QuasiCyclic where

import Data.Bit
import Data.Bits
import Data.Matrix
import Data.Char
import qualified Data.Matrix.Matlab as Matlab

import Data.Word

-- | Not abstract
--
-- You define a n-QuasiCyclic matrix by only including every n-th row.
--
data QuasiCyclic a = QuasiCyclic Int (Matrix a)

toBitMatrix :: Bits a => QuasiCyclic a -> Matrix Bit
toBitMatrix (QuasiCyclic sz a) = fromLists
        [ [ mkBit $ (a ! (m,n)) `testBit` ((j - i) `mod` sz)
          | n <- [1..ncols a], j <- [0..sz-1]
          ]
        | m <- [1..nrows a], i <- [0..sz-1]
        ]

fromBitMatrix :: Bits a => Int -> Matrix Bit -> QuasiCyclic a
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
                    , getBit v  -- if set
                    ]
                  | n <- [1,(sz+1)..ncols a]
                  ]
                | m <- [1,(sz+1)..nrows a]
                ]
          zero = bit 0 `xor` bit 0 -- Hack to get a zero, without


cycleSize :: Bits a => QuasiCyclic a -> Int
cycleSize (QuasiCyclic n a) = n

instance Show a => Show (QuasiCyclic a) where
        -- The first line is the cycle size; the rest is a Matlab style matrix
        show (QuasiCyclic sz a) = show sz ++ "\n" ++ show (Matlab.Matlab a)

instance Read a => Read (QuasiCyclic a) where
        readsPrec _ txt0 = [ (QuasiCyclic sz a,txt2)
                           | (sz,txt1)              <- reads (dropWhile isSpace txt0)
                           , (Matlab.Matlab a,txt2) <- reads (dropWhile isSpace txt1)
                           ]


example :: QuasiCyclic Word8
example = QuasiCyclic 8 $ fromLists [[ 0x01, 0xf8], [0x0, 0xd1]]

bitM :: Matrix Bit
bitM = fromLists [[0,1,0,1,1,0],[0,0,1,0,1,1],[1,0,0,1,0,1]]
