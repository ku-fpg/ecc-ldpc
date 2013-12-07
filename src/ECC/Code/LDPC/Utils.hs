{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Utils where

import Data.Bit
import ECC.Types
import Data.BitMatrix.Loader

mkLDPC :: (MatrixLoader g, MatrixLoader h)
       => String -> String -> Int
       -> (g -> [Bit] -> IO [Bit])
       -> (h -> Int -> [Double] -> IO [Bit])
       -> IO ECC
mkLDPC prefix codeName maxI encoder decoder = do
   g <- loadMatrix (codeName ++ "/G")   -- with G, we prepend the identity
   print (getNRows g, getNCols g)
   h <- loadMatrix (codeName ++ "/H")
   print (getNRows h, getNCols h)
   return $ ECC
        { name     = "ldpc/" ++ prefix ++ "/" ++ codeName ++ "/" ++ show maxI
        , encode   = encoder g
        , decode   = \ inp -> do
                             res <- decoder h maxI inp
                             return $ (,True) $ take (getNRows g) $ res
        , message_length  = getNRows g
        , codeword_length =  getNRows g + getNCols g
        }

-- Our version of atanh
atanh' :: Double -> Double
atanh' x | isInfinite y = signum x * 18.714973875118524
         | otherwise    = y
 where y = atanh x
