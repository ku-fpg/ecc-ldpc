{-# LANGUAGE GADTs, DataKinds, KindSignatures, StandaloneDeriving, TypeFamilies, FlexibleInstances #-}
module Data.BitMatrix.Sparse where

import Data.Bit
import qualified Data.BitVector.Sparse as BV
import Data.BitVector.Sparse (BitVector)
import qualified Data.Vector as V
import Data.Vector (Vector)
import qualified Data.List as L
import Data.BitMatrix.Alist
import Control.Monad.State.Lazy
{-
import Data.Array
import Data.Bit
import Data.Bits
import Data.List
import Data.Alist

import Codes.MoonLDPC

import System.Random
import Debug.Trace
import Data.Array.Matrix (S(..))

-}

data C = CRS | CCS      -- the compression (row or column)

type family T (a :: C) :: C

type instance T CRS = CCS
type instance T CCS = CRS

-- Simple bit-matrix. Use a phantom for CRS vs CCS
data BitMatrix :: C -> * where
  BitMatrixCRS :: Int   -- rows
               -> Int   -- cols
               -> Vector BitVector
               -> BitMatrix CRS
  BitMatrixCCS :: Int   -- rows
               -> Int   -- cols
               -> Vector BitVector
               -> BitMatrix CCS

-- deriving instance Show (BitMatrix a)


class Compress c where
  coerceCompress :: BitMatrix c' -> BitMatrix c

instance Compress CRS where
  coerceCompress = crs

instance Compress CCS where
  coerceCompress = ccs

fromLists :: [[Bit]] -> BitMatrix CRS
fromLists xss = BitMatrixCRS n m
              $ V.fromList $ map BV.fromList xss
  where
     (n,m) = (length xss, length (head xss))

instance Show (BitMatrix a) where
  show (BitMatrixCRS n m v) = unlines [ take m (show (v V.! (i-1)) ++ repeat '0') | i <- [1..n]]
  show bm = unlines $ L.transpose $ lines $ show $ transpose bm

transpose :: BitMatrix c -> BitMatrix (T c)
transpose (BitMatrixCRS n m r) = BitMatrixCCS m n r
transpose (BitMatrixCCS n m r) = BitMatrixCRS m n r

nrows :: BitMatrix a -> Int
nrows (BitMatrixCRS n m r) = n
nrows (BitMatrixCCS n m r) = n

ncols :: BitMatrix a -> Int
ncols (BitMatrixCRS n m r) = m
ncols (BitMatrixCCS n m r) = m

multStd :: BitMatrix CRS -> BitMatrix CCS -> BitMatrix CRS
multStd (BitMatrixCRS n k1 a1) (BitMatrixCCS k2 m a2) | k1 == k2 = fromLists xss
    where
         xss = [ [  parity $ BV.zipWith (*) x1 x2
                 | x1 <- V.toList a1
                 ]
              | x2 <- V.toList a2
              ]
         parity bv = if odd $ length $ BV.elems bv then 1 else 0

ccs :: BitMatrix c -> BitMatrix CCS
ccs b@(BitMatrixCCS {})   = b
ccs (BitMatrixCRS n m arr) = BitMatrixCCS n m
                           $ fmap BV.vector
                           $ V.accum (flip (:)) (V.replicate m [])
                           $ [ (v-1,i+1)
                             | (vs,i) <- V.toList arr `zip` [0..]
                             , v <- BV.elems vs
                             ]

crs :: BitMatrix c -> BitMatrix CRS
crs b@(BitMatrixCRS {})   = b
crs (BitMatrixCCS n m arr) = BitMatrixCRS n m
                           $ fmap BV.vector
                           $ V.accum (flip (:)) (V.replicate n [])
                           $ [ (v-1,i+1)
                             | (vs,i) <- V.toList arr `zip` [0..]
                             , v <- BV.elems vs
                             ]

{-
instance Compress a => Alist (BitMatrix a) where
--  readAlist :: String -> IO (BitMatrix 'CRS)
  readAlist fileName = do txt <- readFile fileName
                          return $ evalState parser $ filter (/= 0) $ map read $ words txt
    where
      parser = do
        n <- item
        m <- item
        _ <- item
        _ <- item
        num_nlist <- replicateM n item
        num_mlist <- replicateM m item
        nlists <- sequence [ replicateM c item | c <- num_nlist ]
        mlists <- sequence [ replicateM c item | c <- num_mlist ]
--        let m2 = rowBitMatrix n (map BitVector mlists)
        let m1 = BitMatrixCRS n m (V.fromList (map BV.vector nlists))
        return (coerceCompress m1)
  -- TODO: add the zeros

--  writeAlist :: String -> (BitMatrix a) -> IO ()
  writeAlist fileName m = writeFile fileName $ unlines $ map unwords $ map (map show) $
        [ [ nrows m, ncols m ]
        , [ maximum $ num_nlist, maximum $ num_mlist ]
        , num_nlist
        , num_mlist
        ]  ++ map BV.elems (V.toList m_crs)
           ++ map BV.elems (V.toList m_ccs)
     where
          BitMatrixCRS _ _  m_crs = crs m
          BitMatrixCCS _ _  m_ccs = ccs m

          num_nlist = map length $ map BV.elems $ V.toList m_crs
          num_mlist = map length $ map BV.elems $ V.toList m_ccs

item :: State [Int] Int
item = do (x:xs) <- get
          put xs
          return x

-}
