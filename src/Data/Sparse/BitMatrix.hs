module Data.Sparse.BitMatrix
  (SparseBitM
  ,toSparseBitM
  ,sparseSet
  ,isSet
  ,sparseBitNrows
  ,sparseBitNcols
  ,(*|)
  )
  where

import qualified Data.Set as S
import qualified Data.Vector.Unboxed as U
import qualified Data.Matrix as M

import Data.Maybe (catMaybes)

type V a = U.Vector a
type M a = M.Matrix a

data SparseBitM =
  SparseBitM
    !Int -- # rows
    !Int -- # cols
    (S.Set (Int, Int))

toSparseBitM :: M Bool -> SparseBitM
toSparseBitM mat
  = SparseBitM (M.nrows mat) (M.ncols mat)
  . S.fromList
  . catMaybes
  $ (go <$> [1..M.nrows mat]
        <*> [1..M.ncols mat])
  where
    go m n
      | mat M.! (m, n) == False = Nothing
      | otherwise               = Just (m, n)

sparseSet :: SparseBitM -> S.Set (Int, Int)
sparseSet (SparseBitM _ _ s) = s

isSet :: SparseBitM -> (Int, Int) -> Bool
isSet sp loc = loc `S.member` sparseSet sp

sparseBitNrows, sparseBitNcols :: SparseBitM -> Int
sparseBitNrows (SparseBitM r _ _) = r
sparseBitNcols (SparseBitM _ c _) = c

(*|) :: (Eq a, Show a, Num a, U.Unbox a) =>
  SparseBitM -> V a -> V Bool
(*|) mat vec = U.fromList (map (fromBit . go) [1..sparseBitNrows mat])
  where
    go r
      = U.sum
      $ U.ifilter (\c _ -> (r, c) `S.member` (sparseSet mat))
                  vec
    {-# INLINE go #-}
{-# INLINE (*|) #-}

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

