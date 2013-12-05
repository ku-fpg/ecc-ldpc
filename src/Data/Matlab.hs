{-# LANGUAGE FlexibleInstances #-}
--
-- | We can read and write matrixes in "Matlab" format,
-- which is just rows of numbers, seperated by spaces,
-- and each row is terminated by a newline.

module Data.Matlab where

import Data.Matrix as M
import Data.Bit

class Matlab a where
  readMatlab  :: String -> IO a
  writeMatlab :: String -> a -> IO ()

instance Matlab (Matrix Bit) where
  readMatlab fileName = do
          txt <- readFile fileName
          return $ M.fromLists $ map (map read) $ map words $ lines txt
