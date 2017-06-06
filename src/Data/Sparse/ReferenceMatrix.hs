-- This supports all the operations of the other ones, but it uses
-- a typical matrix representation.

{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE ScopedTypeVariables #-}

module Data.Sparse.ReferenceMatrix
  (RefMatrix
  ,toSparseBitM
  ,sparseSet
  ,Data.Sparse.ReferenceMatrix.getCol
  ,Data.Sparse.ReferenceMatrix.fromList
  ,(*|)
  ,(!!!)
  ,isSet
  ,Data.Sparse.ReferenceMatrix.nrows
  ,Data.Sparse.ReferenceMatrix.ncols
  ,zero
  )
  where

import qualified Data.Set as S
import Data.Set (Set)
import qualified Data.Matrix as M
import Data.Matrix (Matrix)
import qualified Data.Vector.Unboxed as U

import Data.Coerce
import Data.List (lookup)

import Data.Bit
import Data.Maybe (catMaybes)

type V a = U.Vector a

newtype RefMatrix a = RefMatrix (Matrix a)
  deriving (Functor)

toSparseBitM :: Matrix Bool -> RefMatrix Bit
toSparseBitM = fmap toBit . coerce

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
      $ U.ifilter (\c _ -> mat !!! (r, c+1) == 1)
                  vec
    {-# INLINE go #-}
{-# INLINE (*|) #-}

(!!!) :: (Num a, Eq a) => RefMatrix a -> (Int, Int) -> a
m !!! coords = coerce m M.! coords

isSet :: RefMatrix Bit -> (Int, Int) -> Bool
isSet m coords = m !!! coords == 1

sparseSet :: RefMatrix Bit -> Set (Int, Int)
sparseSet m = S.fromList (catMaybes (go <$> [1..nrows m] <*> [1..ncols m]))
  where
    go r c
      | m !!! (r, c) == 1 = Just (r, c)
      | otherwise         = Nothing

nrows, ncols :: RefMatrix a -> Int
nrows = M.nrows . (coerce :: RefMatrix b -> Matrix b)
ncols = M.ncols . (coerce :: RefMatrix b -> Matrix b)

zero :: forall a. Num a => Int -> Int -> RefMatrix a
zero r c = coerce (M.zero r c :: Matrix a)

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

toBit :: Num a => Bool -> a
toBit False = 0
toBit True  = 1

