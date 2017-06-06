-- This supports all the operations of the other ones, but it uses
-- a typical matrix representation.

module Data.Sparse.ReferenceMatrix
  (RefMatrix
  ,Data.Sparse.ReferenceMatrix.getCol
  ,Data.Sparse.ReferenceMatrix.fromList
  ,(*|)
  ,(!!!)
  ,isSet
  ,Data.Sparse.ReferenceMatrix.nrows
  ,Data.Sparse.ReferenceMatrix.ncols
  )
  where

import qualified Data.Matrix as M
import Data.Matrix (Matrix)
import qualified Data.Vector.Unboxed as U

import Data.Coerce
import Data.List (lookup)

import Data.Bit

type V a = U.Vector a

newtype RefMatrix a = RefMatrix (Matrix a)

getCol :: U.Unbox a => Int -> RefMatrix a -> V a
getCol c = U.convert . M.getCol c . coerce

fromList :: Num a => Int -> Int -> [((Int, Int), a)] -> RefMatrix a
fromList r c a = coerce $ M.matrix r c (`defLookup` a)
  where
    defLookup x y =
      case lookup x y of
        Just r -> r
        _      -> 0

(*|) :: (Eq a, Show a, Num a, U.Unbox a) =>
  RefMatrix Bit -> V a -> V Bool
(*|) mat vec = U.fromList (map (fromBit . go) [1..nrows mat])
  where
    go r
      = U.sum
      $ U.ifilter (\c _ -> mat !!! (r, c) == 1)
                  vec
    {-# INLINE go #-}
{-# INLINE (*|) #-}

(!!!) :: (Num a, Eq a) => RefMatrix a -> (Int, Int) -> a
m !!! coords = coerce m M.! coords

isSet :: RefMatrix Bit -> (Int, Int) -> Bool
isSet m coords = m !!! coords == 1


nrows, ncols :: RefMatrix a -> Int
nrows = M.nrows . (coerce :: RefMatrix b -> Matrix b)
ncols = M.ncols . (coerce :: RefMatrix b -> Matrix b)

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

