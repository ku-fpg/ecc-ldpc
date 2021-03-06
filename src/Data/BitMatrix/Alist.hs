{-# LANGUAGE FlexibleInstances #-}
-- | An Alist is a standard way of defining sparse bit matrixes.
--
-- This is the format defined in
-- <http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html>
--
module Data.BitMatrix.Alist (Alist(..)) where

import Control.Monad.State.Lazy as S
import qualified Data.Vector.Unboxed as U
import Data.Matrix

newtype Alist = Alist { unAlist :: Matrix Bool }

instance Show Alist where
    show (Alist mx) = unlines $ map unwords $ map (map show) $
        [ [ nrows $ mx, ncols $ mx ]
        , [ maximum $ num_nlist, maximum $ num_mlist ]
        , num_nlist
        , num_mlist
        ]  ++ m_crs
           ++ m_ccs
     where
          m_crs = [ [ m | m <- [1..ncols mx], mx!(n,m) ] | n <- [1..nrows mx] ]
          m_ccs = [ [ n | n <- [1..nrows mx], mx!(n,m) ] | m <- [1..ncols mx] ]

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
        return $ fromList n m
               $ [ i `elem` ns
                 | ns <- nlists -- all the rows
                 , i <- [1..m]
                 ]

  -- TODO: add the zeros

item :: State [Int] Int
item = do (x:xs) <- S.get
          S.put xs
          return x

