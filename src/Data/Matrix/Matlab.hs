--
-- | We can read and write matrixes in "Matlab" format,
-- which is just rows of numbers, seperated by spaces,
-- and each row is terminated by a newline.

module Data.Matrix.Matlab where

import Data.Matrix as M
import Data.Bit

newtype Matlab a = Matlab { unMatlab :: Matrix a }

instance Show a => Show (Matlab a) where
    show (Matlab mx) = unlines
        [ unwords [ show (mx ! (m,n)) | n <- [1..ncols mx] ]
        | m <- [1.. nrows mx]
        ]

instance Read a => Read (Matlab a) where
    readsPrec _ txt = [ (Matlab $ M.fromLists $ map (map read) $ map words $ lines txt,"") ]

