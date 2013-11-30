{-# LANGUAGE GADTs, DataKinds, KindSignatures, StandaloneDeriving, TypeFamilies, FlexibleInstances #-}
module Data.BitMatrix.Sparse where

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
import Control.Monad.State.Lazy
-}

{-
newtype BitVector = BitVector { unBitVector :: [Int] }
        deriving (Eq, Ord, Show)

bitVector :: [Bit] -> BitVector
bitVector vs = BitVector [ i | (i,1) <- [1..] `zip` vs ]

--sshowBV :: BitVector -> String
--showBV (BitVector bv) = ""
showBV (BitVector bv) = concat [ take (n-1) (repeat '0') ++ "1" | n <- zipWith (-) bv (0 : bv) ]

zipWithBitVector :: (Bit -> Bit -> Bit) -> BitVector -> BitVector -> BitVector
zipWithBitVector f _ _ | f 0 0 /= 0 = error "zipWithBitVector assumes f 0 0 == 0"
zipWithBitVector f (BitVector xs) (BitVector ys) = BitVector [ z | (z,1) <- loop xs ys ]
  where
        loop []     []     = []
        loop (x:xs) []     = (x,f 1 0) : loop xs []
        loop []     (y:ys) = (y,f 0 1) : loop ys []
        loop (x:xs) (y:ys)
                | x == y   = (x,f 1 1) : loop xs ys
                | x < y    = (x,f 1 0) : loop xs (y:ys)
                | x > y    = (y,f 0 1) : loop (x:xs) ys

data C = CRS | CCS

type family T (a :: C) :: C

type instance T CRS = CCS
type instance T CCS = CRS

-- Simple bit-matrix. Use a phantom for CRS vs CCS
data BitMatrix :: C -> * where
  BitMatrixCRS :: (Int,Int)  -- rows x cols
               -> Array Int BitVector
               -> BitMatrix CRS
  BitMatrixCCS :: (Int,Int)  -- rows x cols
               -> Array Int BitVector
               -> BitMatrix CCS

deriving instance Show (BitMatrix a)

class Compress c where
  coerceCompress :: BitMatrix c' -> BitMatrix c

instance Compress CRS where
  coerceCompress = crs

instance Compress CCS where
  coerceCompress = ccs

showBM :: BitMatrix a -> String
showBM (BitMatrixCRS (n,m) a) =
        unlines [ take m (showBV (a ! i) ++ repeat '0') | i <- [1..n]]
showBM bm = unlines $ transpose $ lines $ showBM $ transposeBM bm

rowBM :: Int -> BitVector ->BitMatrix CRS
rowBM s bv = BitMatrixCRS
                (1,s)
                (array (1,1) [(1,bv)])

columnBM :: Int -> BitVector -> BitMatrix CCS
columnBM s = transposeBM . rowBM s

bitMatrixCRS :: [[Bit]] -> BitMatrix CRS
bitMatrixCRS xss = BitMatrixCRS
                (n,m)
                (array (1,n)
                  [ (i, bitVector xs)
                  | (i,xs) <- [1..] `zip` xss
                  ])
  where
     (n,m) = (length xss, length (head xss))

rowBitMatrix :: Int -> [BitVector] -> BitMatrix CRS
rowBitMatrix w vs = BitMatrixCRS
                (length vs,w)
                (array (1,length vs)
                        [ (i,v)
                        | (i,v) <- [1..] `zip` vs
                        ])


bitMatrixCCS :: [[Bit]] -> BitMatrix CCS
bitMatrixCCS = transposeBM . bitMatrixCRS . transpose


transposeBM :: BitMatrix c -> BitMatrix (T c)
transposeBM (BitMatrixCRS (n,m) r) = BitMatrixCCS (m,n) r
transposeBM (BitMatrixCCS (n,m) r) = BitMatrixCRS (m,n) r

sizeBM :: BitMatrix c -> (Int,Int)
sizeBM (BitMatrixCRS s _) = s
sizeBM (BitMatrixCCS s _) = s

-- force ccs
ccs :: BitMatrix c -> BitMatrix CCS
ccs b@(BitMatrixCCS {})   = b
ccs (BitMatrixCRS (n,m) arr) =
        BitMatrixCCS
                (n,m)
                (fmap (BitVector . reverse)
                  $ accumArray (flip (:)) [] (1,m)
                  [ (j,i)
                  | i <- [1..n]
                  , j <- case arr ! i of BitVector vs -> vs
                  ])
crs :: BitMatrix c -> BitMatrix CRS
crs b@(BitMatrixCRS {})   = b
crs (BitMatrixCCS (n,m) arr) =
        BitMatrixCRS
                (n,m)
                (fmap (BitVector . reverse)
                  $ accumArray (flip (:)) [] (1,n)
                  [ (j,i)
                  | i <- [1..m]
                  , j <- case arr ! i of BitVector vs -> vs
                  ])



mm :: BitMatrix CRS -> BitMatrix CCS -> BitMatrix CRS
mm (BitMatrixCRS (n,k1) a1) (BitMatrixCCS (k2,m) a2) | k1 == k2 = bitMatrixCRS xss
    where
         xss = [ [  parity $ zipWithBitVector (*) (a1 ! i) (a2 ! j)
                 | j <- [1..m]
                 ]
              | i <- [1..n]
              ]
         parity (BitVector xs) = if odd $ length xs then 1 else 0

m1 = bitMatrixCRS [[1,0]]
m1' = columnBM 3 (bitVector [1,0,1])
m2 = bitMatrixCCS [[1,0,1],[0,1,1]]
m2' = bitMatrixCRS [[1,0,1],[0,1,1]]

--constant :: BitMatrix CRS -> (Int,Int) -> Int -> BitVector -> Maybe [Int]
constant m (x,y) k vs | trace (show ("constant",(x,y),k,vs) ++ "\n" ++ showBM m) False = undefined
constant (BitMatrixCRS (n,m) a) (x,y) w bv@(BitVector vs)
  | matches == BitVector [] = Nothing        -- all done
  | otherwise = Just $
        [ (i,matches)
        | i <- [1..n] -- [x+1..n]
        , x /= i               -- XOR by yourself is bad
        , chop y (y+w-1) (a ! i) == matches
        , (a ! i) /= (a ! x)
        ]
  where
        -- what we are looking for
        matches :: BitVector
        matches = zipWithBitVector xor (chop y (y+w-1) (a ! x)) bv

-- xor the 1st and 2nd argument rows, putting the answer in the first argument
xorRow :: Int -> Int -> BitMatrix CRS -> BitMatrix CRS
xorRow r1 r2 (BitMatrixCRS (n,m) a) | trace (show ("xor",r1,r2,S $ showBV $ a ! r1,S $ showBV $ a ! r2)) False = undefined
xorRow r1 r2 (BitMatrixCRS (n,m) a) = sanity $ BitMatrixCRS (n,m)
                                    (a // [ (r1,zipWithBitVector xor (a ! r1) (a ! r2)) ])


sanity = id -- bm@(BitMatrixCRS (n,_) a) | n == (length $ group $ sort $ elems $ a) = bm


findG s h = fn h xs
  where
     (n,m) = sizeBM h

     h' = fn

     fn m ((x,y):xys) = trace (show (x,y)) $
                       case constant m (x,1 + s) y (if x <= y then BitVector [x] else BitVector []) of
                         Nothing -> fn m xys
                         Just [] -> [m] -- error $ show (x,y) ++ "\n" ++ showBM m
                         Just rs -> [ r | (x',_) <- rs, r <- fn (xorRow x x' m) xys ]
     fn m _ = [m]
     xs =
          [ (x,y)
          | y <- [1..n]
          , x <- [1..n]
          , x >= y
          ] ++
          [ (x,y)
          | x <- [1..n]
          , y <- [1..n]
          , x < y
          ]

--     s = 60
--     s = m - n


chop :: Int -> Int -> BitVector -> BitVector
chop lo hi (BitVector xs) = BitVector $ map (\ x -> 1 + x - lo) $ takeWhile (<= hi) $ dropWhile (< lo) $ xs
--stripe m (x,y) = [ m ! (x,y) | y <- [y..columns m]]

rndMatrix = bitMatrixCRS $ unfoldr f xs
  where
        (n,m) = (40,80)
        f [] = Nothing
        f xs = Just (take m xs,drop m xs)
        g = mkStdGen 100
        xs = map (\ x -> if x < 0.09 then 1 else 0) $ take (n * m) $ (randoms g :: [Float])

test = constant rndMatrix (4,11) 3 (BitVector [1,3])

m = bitMatrixCRS
  $ fmap (fmap (\ c -> if c == '0' then 0 else 1))
--  $ reverse
    ["010010",
     "100010",
     "001001",
     "110101",
     "100111",
     "100011"]

--------------------------------------------------------------

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
        let m1 = rowBitMatrix m (map BitVector nlists)
        return (coerceCompress m1)
  -- TODO: add the zeros
--  writeAlist :: String -> (BitMatrix a) -> IO ()
  writeAlist fileName m = writeFile fileName $ unlines $ map unwords $ map (map show) $
        [ [ fst $ nm, snd $ nm ]
        , [ maximum $ num_nlist, maximum $ num_mlist ]
        , num_nlist
        , num_mlist
        ]  ++ [ vs | BitVector vs <- elems m_crs ]
           ++ [ vs | BitVector vs <- elems m_ccs ]
     where
          BitMatrixCRS nm  m_crs = crs m
          BitMatrixCCS nm' m_ccs = ccs m

          num_nlist = map (\ (BitVector vs) -> length vs) $ elems m_crs
          num_mlist = map (\ (BitVector vs) -> length vs) $ elems m_ccs




item :: State [Int] Int
item = do (x:xs) <- get
          put xs
          return x

-- to call (!) later, also generalize to any CCS/CRS
idx :: BitMatrix CRS -> (Int,Int) -> Bit
idx (BitMatrixCRS (n,m) a) (x,y) = if inRange ((1,1),(n,m)) (x,y)
              then case dropWhile (< y) (unBitVector (a ! x)) of
                     (r:_) | r == y -> 1
                     _              -> 0
              else 0
-}