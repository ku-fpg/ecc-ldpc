{-# LANGUAGE FlexibleInstances #-}
-- | An Alist is a standard way of defining sparse bit matrixes.
--
-- This is the format defined in
-- <http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html>
--
module Data.BitMatrix.Alist (Alist(..)) where

import Control.Monad.State.Lazy as S
import Data.Matrix.Unboxed hiding (map, sequence)
import qualified Data.Vector.Unboxed as U
import Data.Bit

newtype Alist = Alist { unAlist :: Matrix Bool }

instance Show Alist where
    show (Alist mx) = unlines $ map unwords $ map (map show) $
        [ [ rows $ mx, cols $ mx ]
        , [ maximum $ num_nlist, maximum $ num_mlist ]
        , num_nlist
        , num_mlist
        ]  ++ m_crs
           ++ m_ccs
     where
          m_crs = [ [ m | m <- [1..cols mx], mx!(n,m) ] | n <- [1..rows mx] ]
          m_ccs = [ [ n | n <- [1..rows mx], mx!(n,m) ] | m <- [1..cols mx] ]

          num_nlist = map length m_crs
          num_mlist = map length m_ccs

instance Read Alist where
  readsPrec _ txt = (\ a -> [(a,"")])  . Alist . evalState parser . filter (/= 0) . filter (/= 0) $ map read $ words txt
    where
      parser = do
        n <- item
        m <- item
        _ <- item
        _ <- item
        num_nlist <- replicateM n item
        num_mlist <- replicateM m item
        nlists <- sequence [ replicateM c item | c <- num_nlist ]
        _mlists <- sequence [ replicateM c item | c <- num_mlist ]
        return $ fromVector (n, m)
               $ U.fromList
               $ [ i `elem` ns
                 | ns <- nlists -- all the rows
                 , i <- [1..m]
                 ]

  -- TODO: add the zeros

item :: State [Int] Int
item = do (x:xs) <- S.get
          S.put xs
          return x

