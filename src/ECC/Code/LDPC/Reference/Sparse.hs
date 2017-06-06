{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TupleSections #-}

-- Uses a sparse bit matrix and a sparse matrix for soft values

module ECC.Code.LDPC.Reference.Sparse where

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix as M
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V

import Data.Map (Map)
import qualified Data.Map as Map

import Data.Set (Set)
import qualified Data.Set as Set

import Data.Maybe (catMaybes)
import Control.Applicative hiding ((<|>))
import Data.Foldable as F

import Data.List (groupBy)
import Data.Function (on)

type M a = Matrix a
type V a = U.Vector a

type SparseM a  = Map (Int, Int) a
data SparseBitM =
  SparseBitM
    !Int -- # rows
    !Int -- # cols
    (Set (Int, Int))

sparseGetCol :: U.Unbox a => Int -> SparseM a -> U.Vector a
sparseGetCol c sp = U.fromList . Map.elems $ Map.filterWithKey (\(_, c') _ -> c' == c) sp
{-# INLINE sparseGetCol #-}

toSparseM :: (Eq a, Num a) => M a -> SparseM a
toSparseM mat =
    Map.fromList $ catMaybes (go <$> [1..nrows mat] <*> [1..ncols mat])
  where
    go r c
      | v == 0    = Nothing
      | otherwise = Just ((r, c), v)
      where
        v = mat ! (r, c)

toSparseBitM :: M Bool -> SparseBitM
toSparseBitM mat
  = SparseBitM (M.nrows mat) (M.ncols mat)
  . Set.fromList
  . catMaybes
  $ (go <$> [1..M.nrows mat]
        <*> [1..M.ncols mat])
  where
    go m n
      | mat ! (m, n) == False = Nothing
      | otherwise             = Just (m, n)

sparseSet :: SparseBitM -> Set (Int, Int)
sparseSet (SparseBitM _ _ s) = s

isSet :: SparseBitM -> (Int, Int) -> Bool
isSet sp loc = loc `Set.member` sparseSet sp

sparseNrows, sparseNcols :: SparseBitM -> Int
sparseNrows (SparseBitM r _ _) = r
sparseNcols (SparseBitM _ c _) = c

(*|) :: (Eq a, Show a, Num a, U.Unbox a) =>
  SparseBitM -> V a -> V Bool
(*|) mat vec = U.fromList (map (fromBit . go) [1..sparseNrows mat])
  where
    go r
      = U.sum
      $ U.ifilter (\c _ -> (r, c) `Set.member` (sparseSet mat))
                  vec
    {-# INLINE go #-}
{-# INLINE (*|) #-}

-- If not found, return zero
(!!!) :: (Num a, Eq a) => SparseM a -> (Int, Int) -> a
m !!! (r, c) =
  case Map.lookup (r, c) m of
    Just v -> v
    _      -> 0

code :: Code
code = mkLDPC_Code "reference-sparse" encoder decoder

---------------------------------------------------------------------

encoder :: M Bool -> U.Vector Bool -> IO (U.Vector Bool)
encoder g v = return $ U.map toBool $ U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) (((identity (nrows g)) <|> fmap fromBool g))))

---------------------------------------------------------------------

decoder :: M Bool -> Int -> U.Vector Double -> IO (U.Vector Bool)
decoder a maxIterations orig_lam = return $ U.convert (ldpc a maxIterations (U.convert orig_lam))

ldpc :: M Bool -> Int -> V Double -> V Bool
ldpc a0 maxIterations orig_lam = U.map hard $ loop 0 orig_ne orig_lam
  where
    a :: SparseBitM
    a = toSparseBitM a0

    SparseBitM _ _ aSet0 = a

    aSet :: [[(Int, Int)]]
    aSet = groupBy ((==) `on` fst) $ Set.toAscList aSet0

    aList :: [((Int, Int), [Int])]
    aList = concatMap
      (\list@((m,_):_) ->
          let ones = rowOnes m
          in
          map (, ones) list)
      aSet

    orig_ne :: SparseM Double
    orig_ne = Map.empty

    rowOnes :: Int -> [Int]
    rowOnes m =
      [ j
      | j <- [1 .. U.length orig_lam]
      , isSet a (m, j)
      ]

    loop :: Int -> SparseM Double -> V Double -> V Double
    loop !n ne lam
      | U.all (== False) ans = lam
      | n >= maxIterations   = orig_lam
      | otherwise            = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = U.map (toBit . hard) lam

        ans :: V Bool
        ans = a *| c_hat

        ne' :: SparseM Double
        ne' = Map.fromList $ map go aList
          where
            go (coords@(m, n), ones) =
              (coords
              ,-2 * atanh' (product
                   [ tanh (- ((lam U.! (j-1) - ne !!! (m,j)) / 2))
                   | j <- ones
                   , j /= n
                   ])
              )

        lam' :: V Double
        lam' = U.fromList [ (orig_lam U.! (j - 1)) + U.sum (sparseGetCol j ne')
                          | j <- [1 .. U.length lam]
                          ]



-- toBits :: M Bool -> BitMatrix CRS
-- toBits = S.fromLists . map (map toBit) . M.toLists

toBit :: Num a => Bool -> a
toBit False = 0
toBit True  = 1

-- fromBits :: BitMatrix CRS -> M Bool
-- fromBits = M.fromLists . map (map fromBit) . S.toLists

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

