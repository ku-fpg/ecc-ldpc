module Data.Sparse.Matrix
  (SparseM
  ,Data.Sparse.Matrix.getCol
  ,toSparseM
  ,Data.Sparse.Matrix.fromList
  ,(!!!)
  )
  where

import qualified Data.Map as Map
import Data.Map (Map)
import Data.Matrix as M
import qualified Data.Vector.Unboxed as U

import Data.Maybe (catMaybes)

import Data.Coerce

type V a = U.Vector a
type M a = M.Matrix a

newtype SparseM a = SparseM (Map (Int, Int) a)

getCol :: U.Unbox a => Int -> SparseM a -> U.Vector a
getCol c = {-# SCC "sparseGetCol" #-}
  U.fromList . Map.elems . Map.filterWithKey (\(_, c') _ -> c' == c) . coerce
{-# INLINE getCol #-}

toSparseM :: (Eq a, Num a) => M a -> SparseM a
toSparseM mat =
    coerce . Map.fromList $ catMaybes (go <$> [1..nrows mat] <*> [1..ncols mat])
  where
    go r c
      | v == 0    = Nothing
      | otherwise = Just ((r, c), v)
      where
        v = mat ! (r, c)

fromList :: [((Int, Int), a)] -> SparseM a
fromList = coerce . Map.fromList

-- If not found, return zero
(!!!) :: (Num a, Eq a) => SparseM a -> (Int, Int) -> a
m !!! (r, c) =
  case Map.lookup (r, c) (coerce m) of
    Just v -> v
    _      -> 0

