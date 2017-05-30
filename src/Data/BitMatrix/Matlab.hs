{-# LANGUAGE FlexibleInstances #-}
--
-- | We can read and write matrixes in "Matlab" format,
-- which is just rows of numbers, seperated by spaces,
-- and each row is terminated by a newline.

module Data.BitMatrix.Matlab where

import Data.Matrix as M

newtype Matlab = Matlab { unMatlab :: Matrix Bool }

instance Show Matlab where
    show (Matlab mx) = unlines
        [ unwords [ showBit (mx ! (m,n)) | n <- [1..ncols mx] ]
        | m <- [1.. nrows mx]
        ]

instance Read Matlab where
    readsPrec _ txt = [ (Matlab $ M.fromLists $ map (map readBit) $ map words $ lines txt,"") ]

readBit :: String -> Bool
readBit "0" = False
readBit "1" = True
readBit _   = error "readBit: no parse"

showBit :: Bool -> String
showBit False = "0"
showBit True  = "1"

