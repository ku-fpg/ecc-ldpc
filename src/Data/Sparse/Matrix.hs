-- NOTE: Not actually a sparse representation, just has a fast getCol and
-- indexing operation.

{-# LANGUAGE TupleSections #-}

module Data.Sparse.Matrix
  (SparseM
  ,Data.Sparse.Matrix.colSum
  -- ,toSparseM
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

import Data.Maybe (catMaybes)

import Data.Coerce

import Data.List (lookup)
import Data.Bit

type V a = V.Vector a
type M a = M.Matrix a

-- | Column-major
newtype SparseM a = SparseM (V.Vector [(Int, a)])

-- getCol :: U.Unbox a => Int -> SparseM a -> U.Vector a
-- getCol c = {-# SCC "getCol" #-}
--   U.convert . M.getRow c . coerce
-- {-# INLINE getCol #-}

colSum :: Num a => Int -> SparseM a -> a
colSum c (SparseM m) = sum $ map snd (m V.! (c-1))

-- toSparseM :: (Eq a, Num a) => M a -> SparseM a
-- toSparseM = coerce
-- toSparseM mat = coerce mat
--     coerce . Map.fromList $ catMaybes (go <$> [1..nrows mat] <*> [1..ncols mat])
--   where
--     go r c
--       | v == 0    = Nothing
--       | otherwise = Just ((r, c), v)
--       where
--         v = mat ! (r, c)

fromList :: Num a => Int -> Int -> [((Int, Int), a)] -> SparseM a
fromList numRows numCols assocs =
  SparseM $ V.fromList (map go [1..numCols])
  -- SparseM $ V.fromList (catMaybes (_ (go <$> [1..numCols] <$> [1..numRows])))
  where
    assocMap = Map.fromList assocs
    go c =
      catMaybes $ map (\r -> (r,) <$> (Map.lookup (r, c) assocMap)) [1..numRows]
    -- go c r =
    --   (r,) <$> (lookup (r, c) assoc)

sparseMatrix :: Num a => Int -> Int -> ((Int, Int) -> Maybe a) -> SparseM a
sparseMatrix numRows numCols f = {-# SCC "sparseMatrix" #-}
  SparseM $ V.fromList (map go [1..numCols])
  where
    go c =
      catMaybes $ map (\r -> (r,) <$> f (r, c)) [1..numRows]
    {-# INLINE go #-}
{-# INLINE sparseMatrix #-}

(!!!) :: (Num a, Eq a) => SparseM a -> (Int, Int) -> a
(SparseM m) !!! (r, c) = sparseLookup r (m V.! (c-1))
  where
    sparseLookup x y =
      case lookup x y of
        Just v -> v
        _      -> 0
  -- case Map.lookup (r, c) (coerce m) of
  --   Just v -> v
  --   _      -> 0
{-# INLINE (!!!) #-}

-- (*|) :: (Eq a, Show a, Num a) =>
--   SparseM Bit -> V a -> V Bool
-- (*|) (SparseM m) v = undefined

-- swap :: (a, b) -> (b, a)
-- swap (x, y) = (y, x)

