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

import ECC.Code.LDPC.Reference.Orig (encoder)

import Data.Sparse.Matrix as SM
import Data.Sparse.BitMatrix as SM

-- import Data.Sparse.ReferenceMatrix as SM

type M a = Matrix a
type V a = U.Vector a


code :: Code
code = mkLDPC_Code "reference-sparse" encoder decoder

---------------------------------------------------------------------

decoder :: M Bool -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a _ maxIterations orig_lam = Just $ U.convert (ldpc a maxIterations (U.convert orig_lam))

ldpc :: M Bool -> Int -> V Double -> V Bool
ldpc a0 maxIterations orig_lam = U.map hard $ loop 0 orig_ne orig_lam
  where
    a :: SparseBitM
    a = toSparseBitM a0

    aSet0 :: Set (Int, Int)
    aSet0 = sparseSet a

    aSet :: [[(Int, Int)]]
    aSet = groupBy ((==) `on` fst) $ Set.toAscList aSet0

    aList :: [((Int, Int), [Int])]
    aList = concatMap
      (\list@((m,_):_) ->
          let ones = rowOnes m
          in
          map (, ones) list)
      aSet

    aListMap :: Map (Int, Int) [Int]
    aListMap = Map.fromList aList

    orig_ne :: SparseM Double
    orig_ne = SM.fromList (SM.nrows a) (SM.ncols a) []
    -- orig_ne = SM.zero (SM.nrows a) (SM.ncols a) --SM.fromList []

    rowOnes :: Int -> [Int]
    rowOnes m =
      [ j
      | j <- [1 .. U.length orig_lam]
      , isSet a (m, j)
      ]

    loop :: Int -> SparseM Double -> V Double -> V Double
    loop !n ne lam
      -- | U.all (== False) ans = lam
      | allMultFalse a c_hat = lam
      | n >= maxIterations   = orig_lam
      | otherwise            = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = U.map (toBit . hard) lam

        -- ans :: V Bool
        -- ans = a *| c_hat

        ne' :: SparseM Double
        ne' = SM.sparseMatrix (SM.nrows a) (SM.ncols a) go
          where
            go (m, n)
              | isSet a (m, n) =
                Just $ -2 * atanh' (product
                  [ tanh (- ((lam U.! (j-1) - ne !!! (m,j)) / 2))
                  | j <- ones
                  , j /= n
                  ])
              | otherwise = Nothing
              where
                Just ones = Map.lookup (m, n) aListMap


        -- ne' :: SparseM Double
        -- ne' = SM.fromList (SM.nrows a) (SM.ncols a) $ map go aList
        --   where
        --     go (coords@(m, n), ones) =
        --       (coords
        --       ,-2 * atanh' (product
        --            [ tanh (- ((lam U.! (j-1) - ne !!! (m,j)) / 2))
        --            | j <- ones
        --            , j /= n
        --            ])
        --       )

        lam' :: V Double
        lam' = U.fromList [ (orig_lam U.! (j - 1)) + (SM.colSum j ne')
                          | j <- [1 .. U.length lam]
                          ]

toBit :: (Num a) => Bool -> a
toBit False = 0
toBit True  = 1

