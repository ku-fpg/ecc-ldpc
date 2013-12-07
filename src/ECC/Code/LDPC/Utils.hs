{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Utils where

import Data.Bit
import ECC.Types
import Data.Matrix
import qualified Data.Vector as V
import Data.Vector(Vector)
import Data.Alist
import Data.BitMatrix.Loader


type G = Matrix Bit
type H = Matrix Bit

mkLDPC :: (MatrixLoader g, MatrixLoader h)
       => String -> String -> Int
       -> (g -> [Bit] -> IO [Bit])
       -> (h -> Int -> [Double] -> IO [Bit])
       -> IO ECC
mkLDPC prefix codeName maxI encoder decoder = do
   g <- loadMatrix (codeName ++ "/G")   -- with G, we prepend the identity
   h <- loadMatrix (codeName ++ "/H")
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
