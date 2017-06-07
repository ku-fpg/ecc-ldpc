-- NOTE: Not actually a sparse representation, just has a fast getCol and
-- indexing operation.

{-# LANGUAGE TupleSections #-}

module Data.Sparse.Matrix
  (SparseM
  ,Data.Sparse.Matrix.colSum
  ,Data.Sparse.Matrix.fromList
  ,sparseMatrix
  ,(!!!)
  )
  where

import qualified Data.Map as Map
import Data.Map (Map)
import Data.Matrix as M
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Set as S

import Data.Maybe (catMaybes)

import Data.Coerce

import Data.List (lookup)
import Data.Bit

type V a = V.Vector a
type M a = M.Matrix a

-- | Column-major
newtype SparseM a = SparseM (V.Vector [(Int, a)])

colSum :: Num a => Int -> SparseM a -> a
colSum c (SparseM m) = sum $ map snd (m V.! (c-1))

fromList :: Num a => Int -> Int -> [((Int, Int), a)] -> SparseM a
fromList numRows numCols assocs =
  SparseM $ V.fromList (map go [1..numCols])
  where
    assocMap = Map.fromList assocs
    go c =
      catMaybes $ map (\r -> (r,) <$> (Map.lookup (r, c) assocMap)) [1..numRows]


-- | NOTE: Expects assoc list to be (c, r) format
sparseMatrix :: Num a => Int -> Int -> [(Int, [Int])] -> ((Int, Int) -> a) -> SparseM a
sparseMatrix numRows numCols indices f = {-# SCC "sparseMatrix" #-}
  SparseM $ V.fromList (map go indices)
  where
    go (c, rs) = map (\r -> (r, f (r, c))) rs
{-# INLINE sparseMatrix #-}

(!!!) :: (Num a, Eq a) => SparseM a -> (Int, Int) -> a
(SparseM m) !!! (r, c) = sparseLookup r (m V.! (c-1))
  where
    sparseLookup x y =
      case lookup x y of
        Just v -> v
        _      -> 0
{-# INLINE (!!!) #-}

