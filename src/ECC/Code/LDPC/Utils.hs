{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Utils where

import Data.Bit
import ECC.Types
import Data.Matrix
import qualified Data.Vector as V
import Data.Vector(Vector)
import Data.Alist

type G = Matrix Bit
type H = Matrix Bit

mkLDPC :: String -> String -> Int -> (G -> Vector Bit -> Vector Bit) -> (H -> Int -> Vector Double -> IO (Vector Bit)) -> IO ECC
mkLDPC prefix codeName maxI encoder decoder = do
   g :: G <- readAlist ("codes/" ++ codeName ++ ".G")   -- with G, we prepend the identity
   let full_g = identity (nrows g) <|> g
   h :: H <- readAlist ("codes/" ++ codeName ++ ".H")
   return $ ECC
        { name     = "ldpc/" ++ prefix ++ "/" ++ codeName ++ "/" ++ show maxI
        , encode   = return . V.toList . encoder full_g . V.fromList
        , decode   = \ inp -> do
                             res <- decoder h maxI (V.fromList inp)
                             return $ (,True) $  take (nrows full_g) $ V.toList res
        , message_length  = nrows full_g
        , codeword_length =  ncols full_g
        }

-- Our version of atanh
atanh' :: Double -> Double
atanh' x | isInfinite y = signum x * 18.714973875118524
         | otherwise    = y
 where y = atanh x
