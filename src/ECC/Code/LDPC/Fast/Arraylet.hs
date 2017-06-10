{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Fast.Arraylet where

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Matrix.QuasiCyclic as Q
import Debug.Trace

type M a = Matrix a
type V a = U.Vector a

-- An arraylet is a small array, that acts like a NxN matrix, with one element in each row and column,
-- (a rotated identity matrix, for example). We index explicity start at row/column *0*.
data Arraylet a = Arraylet Int (U.Vector a)
 deriving Show

sizeOfArraylet :: Arraylet a -> Int
sizeOfArraylet (Arraylet n _) = n

--indexByRow :: Arraylet a -> Int -> a
--indexByRow (Arraylet a _) = undefined

-- A matrix of arraylets, the same size
data Matrixlet a = Matrixlet Int (Matrix (Maybe (Arraylet a)))
 deriving Show
 
initMatrixlet :: forall a . a -> Q.QuasiCyclic Integer -> Matrixlet a
initMatrixlet zero (Q.QuasiCyclic sz m) = Matrixlet sz (fmap f m)
    where f :: Integer -> Maybe (Arraylet a)
          f 0 = Nothing
          f n = undefined

-- A bit matrix of arraylets.
data BitMatrixlet = BitMatrixlet Int (Matrix (Maybe Int))

code :: Code
code = mkLDPC_Code "arraylet" encoder decoder

---------------------------------------------------------------------

encoder :: M Bool -> Rate -> U.Vector Bool -> U.Vector Bool
encoder g _ v = U.map toBool $ U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) ((fmap fromBool g))))

---------------------------------------------------------------------

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a rate maxIterations orig_lam = Just $ U.convert (ldpc a maxIterations (U.convert orig_lam))

toBit :: Num a => Bool -> a
toBit False = 0
toBit True  = 1

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

multStd' :: M Bool -> M Bool -> M Bool
multStd' a b = fmap fromBit $ multStd (fmap toBit a) (fmap toBit b)

toBits :: M Bool -> M Bit
toBits = fmap toBit

fromBits :: M Bit -> M Bool
fromBits = fmap fromBit

ldpc :: Q.QuasiCyclic Integer -> Int -> V Double -> V Bool
ldpc a maxIterations orig_lam = traceShow msg $ U.map hard $ loop 0 orig_ne orig_lam
  where

    msg = show a

    a' = toBits (Q.toBitMatrix a)

    orig_ne :: M Double
    orig_ne = fmap (const 0) $ Q.toBitMatrix a

    loop :: Int -> M Double -> V Double -> V Double
    loop !n ne lam
        | U.all (== False) ans     = lam
        | n >= maxIterations       = orig_lam
        | otherwise                = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = U.map (toBit . hard) lam

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bool
        ans = U.convert $ fmap fromBit $ getCol 1 (a' `multStd` colVector (U.convert c_hat))

        -- was bug here: V's start at index 0, not 1
        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a' ! (m,n) == 1
                then
                   -- was bug here: we need to cap atanh's answer
                    -2 * atanh' (product
                        [ tanh (- ((lam U.! (j-1) - ne ! (m,j)) / 2))
                        | j <- [1 .. U.length orig_lam]
                        , j /= n
                        , a' ! (m,j) == 1
                        ])
                else 0

        -- Was bug here: needed to add the orig_lam
        lam' :: V Double
        lam' = U.fromList [ U.foldr (+) (orig_lam U.! (j - 1)) (U.convert (getCol j ne'))
                          | j <- [1 .. U.length lam]
                          ]

