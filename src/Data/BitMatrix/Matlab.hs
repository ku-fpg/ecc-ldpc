{-# LANGUAGE FlexibleInstances #-}
--
-- | We can read and write matrixes in "Matlab" format,
-- which is just rows of numbers, seperated by spaces,
-- and each row is terminated by a newline.

module Data.BitMatrix.Matlab where

import Data.Matrix as M
import Data.Bit

newtype Matlab = Matlab { unMatlab :: Matrix Bit }

instance Show Matlab where
    show (Matlab mx) = unlines
        [ unwords [ show (mx ! (m,n)) | m <- [1..ncols mx] ]
        | n <- [1.. nrows mx]
        ]

instance Read Matlab where
    readsPrec _ txt = [ (Matlab $ M.fromLists $ map (map read) $ map words $ lines txt,"") ]

