--
-- | We can read and write matrixes in "Matlab" format,
-- which is just rows of numbers, seperated by spaces,
-- and each row is terminated by a newline.
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE UndecidableInstances #-}

module Data.Matrix.Matlab where

import Data.Matrix.Unboxed (Matrix, (!), cols, rows)
import qualified Data.Matrix.Unboxed as M
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector.Generic as G
import Data.Bit

newtype Matlab a = Matlab { unMatlab :: Matrix a }

instance (G.Vector U.Vector a, Show a) => Show (Matlab a) where
    show (Matlab mx) = unlines
        [ unwords [ show (mx ! (m,n)) | n <- [1..cols mx] ]
        | m <- [1.. rows mx]
        ]

instance (G.Vector U.Vector a, Read a) => Read (Matlab a) where
    readsPrec _ txt = [ (Matlab $ M.fromLists $ map (map read) $ map words $ lines txt,"") ]

