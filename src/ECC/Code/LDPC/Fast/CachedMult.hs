{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, DeriveFunctor #-}
module ECC.Code.LDPC.Fast.CachedMult where

-- Uses the StableDiv data structure to cache multiplications.

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Data.Bit
import Data.Bits
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import qualified Data.Matrix.QuasiCyclic as Q
import Debug.Trace

import qualified ECC.Code.LDPC.Fast.Encoder as E

import Data.Semigroup
import Data.Foldable (foldl')


data StableDiv
  = StableDiv
       !Double -- | The "worst" value for division: the closest to zero
       !Double -- | The result of the multiplications,
               --   excluding the "worst" value

absMinMax :: Double -> Double -> (Double, Double)
absMinMax x y
  | abs x < abs y = (x, y)
  | otherwise     = (y, x)

instance Semigroup StableDiv where
  StableDiv a b <> StableDiv c d =
    StableDiv minOfMins (b * maxOfMins * d)
    where
      (minOfMins, maxOfMins) = absMinMax a c

lit :: Double -> StableDiv
lit x
  | x >= 1    = StableDiv 1 x
  | otherwise = StableDiv x 1

sdiv :: StableDiv -> Double -> Double
sdiv (StableDiv a b) c
  | a == c    = b
  | otherwise = a * (b/c)



type M a = Matrix a
type V a = U.Vector a

-- An arraylet is a small array, that acts like a NxN matrix, with one element in each row and column,
-- (a rotated identity matrix, for example). We index explicity start at row/column *0*.
-- The Int argument is the rotation. The row index is the index into the array.
data Arraylet a = Arraylet Int (U.Vector a)
 deriving (Show)

mapArraylet :: (U.Unbox a, U.Unbox b) => (a -> b) -> Arraylet a -> Arraylet b
mapArraylet f (Arraylet sz a) = Arraylet sz (U.map f a)

sizeOfArraylet :: Arraylet a -> Int
sizeOfArraylet (Arraylet n _) = n

lookupArraylet :: U.Unbox a => Arraylet a -> (Int,Int) -> Maybe a
lookupArraylet (Arraylet off m) (r,c)
     | r < 0 || c < 0           = error "arraylet lookup out of bounds"
     | r >= len || c >= len     = error "arraylet lookup out of bounds"
     | r == (c - off) `mod` len = Just $ m U.! r
     | otherwise                = Nothing
  where len = U.length m

arrayArraylet :: forall a . U.Unbox a => Int -> Int -> ((Int,Int) -> a) -> Arraylet a
arrayArraylet sz off k = Arraylet off $ U.generate sz $ \ r -> k (r,(r + off) `mod` sz)

-- | An indexed map over arraylets (indexing based on 'arrayArraylet')
imapArraylet :: forall a b . (U.Unbox a, U.Unbox b) => Int -> ((Int,Int) -> a -> b) -> Arraylet a -> Arraylet b
imapArraylet sz k (Arraylet off arr) = Arraylet off $ U.imap (\ r v -> k (r,(r + off) `mod` sz) v) arr

foldRowsArraylet :: U.Unbox a => Arraylet a -> U.Vector a
foldRowsArraylet (Arraylet sz m) = b `mappend` a
    where (a,b) = U.splitAt (U.length m - sz) m

foldRowsArraylet' :: U.Unbox a => ((Int,Int) -> a -> b) -> Arraylet a -> V.Vector b
foldRowsArraylet' k (Arraylet off m) = V.generate sz $ \ ix -> let !ix' = (ix + (sz - off)) `mod` sz in
                  k (ix',(ix' + off) `mod` sz) (m U.! ix')
  where !sz = U.length m

foldColsArraylet :: U.Unbox a => Arraylet a -> U.Vector a
foldColsArraylet (Arraylet sz m) = m

foldColsArraylet' :: U.Unbox a => ((Int,Int) -> a -> b) -> Arraylet a -> V.Vector b
foldColsArraylet' k (Arraylet off m) = V.generate sz $ \ ix -> k (ix,(ix + off) `mod` sz) (m U.! ix)
  where !sz = U.length m

--indexByRow :: Arraylet a -> Int -> a
--indexByRow (Arraylet a _) = undefined

-- A matrix of arraylets, all the same size. 0-indexed.
data Matrixlet a = Matrixlet Int (Matrix (Maybe (Arraylet a)))
 deriving Show

initMatrixlet :: forall a . U.Unbox a => a -> Q.QuasiCyclic Integer -> Matrixlet a
initMatrixlet zero (Q.QuasiCyclic sz m) = Matrixlet sz (fmap f m)
    where f :: Integer -> Maybe (Arraylet a)
          f 0 = Nothing
          f n | popCount n == 1 = Just $ Arraylet (g n) (U.replicate sz zero)
              | otherwise       = error $ "QuasiCyclic matrix has non-powers of two initial value of " ++ show n

          -- The number must be a power of two, because there is only one bit set.
          g :: Integer -> Int
          g x | x `testBit` 0 = 0
              | x == 0        = error "got to zero; should never happen"
              | otherwise = 1 + g (x `shiftR` 1)


matrixMatrixlet :: forall a b . U.Unbox b => Matrixlet a -> ((Int,Int) -> b) -> Matrixlet b
matrixMatrixlet (Matrixlet sz a) k = Matrixlet sz $ matrix (nrows a) (ncols a) k'
  where
      k' :: (Int, Int) -> Maybe (Arraylet b)
      k' (r,c) = case a ! (r,c) of
                   Nothing -> Nothing
                   Just (Arraylet off _) -> Just $ arrayArraylet sz off $ \ (r',c') -> k ((r-1) * sz + r',(c-1) * sz + c')

-- | Indexed map over matrixlets (based on 'matrixMatrixlet')
imapMatrixlet :: forall a b . (U.Unbox a, U.Unbox b) => ((Int,Int) -> a -> b) -> Matrixlet a -> Matrixlet b
imapMatrixlet k (Matrixlet sz a) = Matrixlet sz $ matrix (nrows a) (ncols a) k'
  where
      k' :: (Int, Int) -> Maybe (Arraylet b)
      k' (r,c) =
        fmap (imapArraylet sz $ \ (r',c') v ->
                  let r'' = (r-1) * sz + r'
                      c'' = (c-1) * sz + c'
                  in
                  k (r'',c'') v)
             (a ! (r, c))

mapMatrixlet :: (U.Unbox a, U.Unbox b) => (a -> b) -> Matrixlet a -> Matrixlet b
mapMatrixlet f (Matrixlet sz a) = Matrixlet sz (fmap (fmap (mapArraylet f)) a)

-- Assumes there is always one value on every column
foldColsMatrixlet :: U.Unbox a => (a -> a -> a) -> Matrixlet a -> U.Vector a
foldColsMatrixlet f (Matrixlet n a) = U.concat
    [ foldr1 (U.zipWith f) [ foldColsArraylet x | c <- [1..ncols a], Just x <- [a ! (r,c)]]
    | r <- [1..nrows a]
    ]

foldColsMatrixlet' :: U.Unbox a => ((Int,Int) -> a -> b) -> (b -> b -> b) -> Matrixlet a -> V.Vector b
foldColsMatrixlet' f g (Matrixlet sz a) = V.concat
    [ foldr1 (V.zipWith g) [ foldColsArraylet' (\ (r',c') a -> f ((r-1) * sz + r',(c-1) * sz + c') a) x | c <- [1..ncols a], Just x <- [a ! (r,c)]]
    | r <- [1..nrows a]
    ]

foldRowsMatrixlet :: U.Unbox a => (a -> a -> a) -> Matrixlet a -> U.Vector a
foldRowsMatrixlet f (Matrixlet n a) = U.concat
    [ foldr1 (U.zipWith f) [ foldRowsArraylet x | r <- [1..nrows a], Just x <- [a ! (r,c)]]
    | c <- [1..ncols a]
    ]

foldRowsMatrixlet' :: U.Unbox a => ((Int,Int) -> a -> b) -> (b -> b -> b) -> Matrixlet a -> V.Vector b
foldRowsMatrixlet' f g (Matrixlet sz a) = V.concat
    [ foldr1 (V.zipWith g) [ foldRowsArraylet' (\ (r',c') a -> f ((r-1) * sz + r',(c-1) * sz + c') a) x | r <- [1..nrows a], Just x <- [a ! (r,c)]]
    | c <- [1..ncols a]
    ]

toMatrix :: U.Unbox a => a -> Matrixlet a -> Matrix a
toMatrix def mLet@(Matrixlet sz m) = matrix (sz * nrows m) (sz * ncols m) $ \ (x,y) ->
    case mLet `lookupMatrixlet` (x-1,y-1) of
      Nothing -> def
      Just v -> v

lookupMatrixlet :: U.Unbox a => Matrixlet a -> (Int,Int) -> Maybe a
lookupMatrixlet (Matrixlet sz m) (r,c) =
    case m ! (r1+1,c1+1) of
      Nothing -> Nothing
      Just arr -> lookupArraylet arr (r2,c2)

    where (r1,r2) = r `divMod` sz
          (c1,c2) = c `divMod` sz


-- A bit matrix of arraylets.
data BitMatrixlet = BitMatrixlet Int (Matrix (Maybe Int))

code :: Code
code = mkLDPC_Code "arraylet-cm" E.encoder decoder

---------------------------------------------------------------------

decoder :: Q.QuasiCyclic Integer -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a = \ rate maxIterations orig_lam -> Just $ U.convert (ldpc mLet maxIterations (U.convert orig_lam))
  where mLet = initMatrixlet 0 a

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

ldpc :: Matrixlet Double -> Int -> V Double -> V Bool
ldpc mLet maxIterations orig_lam = {- traceShow msg $ -} U.map hard $ loop 0 mLet orig_lam
  where
    cols :: V.Vector [Int]
    cols = foldColsMatrixlet' (\ (r,c) _ -> [c]) (++) mLet

    loop :: Int -> Matrixlet Double -> V Double -> V Double
    loop !n ne lam
        | U.all (== False) ans     = lam
        | n >= maxIterations       = orig_lam
        | otherwise                = loop (n+1) ne' lam'
      where
        ans :: V Bool
        ans = foldColsMatrixlet (/=) $ matrixMatrixlet mLet $ \ (r,c) -> hard (lam U.! c)

        ne_tanh'mat :: Matrixlet Double
        ne_tanh'mat =
          imapMatrixlet
            (\ (m,n) v -> tanh (- ((lam U.! n - v) / 2)))
            ne

        ne_tanhMulted :: V.Vector StableDiv
        ne_tanhMulted = foldColsMatrixlet' (\_ v -> lit v) (<>) ne_tanh'mat

        ne' :: Matrixlet Double
        ne' =
          imapMatrixlet
            (\ (m,n) v -> -2 * atanh' ((ne_tanhMulted V.! m) `sdiv` v))
            ne_tanh'mat

        lam' :: V Double
        lam' = U.zipWith (+) orig_lam $ foldRowsMatrixlet (+) ne'

        unsafeLookup :: Int -> [(Int, a)] -> a
        unsafeLookup i assocs =
          let Just v = lookup i assocs
          in v

