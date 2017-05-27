{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Utils where

import ECC.Types
import ECC.Puncture
import Data.BitMatrix.Loader
import Data.Char (isDigit)
import Data.Semigroup
import qualified Data.Vector.Unboxed as U

mkLDPC :: (MatrixLoader g, MatrixLoader h)
       => String -> String -> Int
       -> (g -> U.Vector Bool -> IO (U.Vector Bool))
       -> (h -> Int -> U.Vector Double -> IO (U.Vector Bool))
       -> IO (ECC IO)
mkLDPC prefix codeName maxI encoder decoder = do
   g <- loadMatrix (codeName ++ "/G")   -- with G, we prepend the identity
   print (getNRows g, getNCols g)
   h <- loadMatrix (codeName ++ "/H")
   print (getNRows h, getNCols h)
   let encoder' = encoder g
   let decoder' = decoder h maxI
   return $ ECC
        { name     = "ldpc/" ++ prefix ++ "/" ++ codeName ++ "/" ++ show maxI
        , encode   = encoder'
        , decode   = \ inp -> do
                             res <- decoder' inp
                             return $ (,True) $ U.take (getNRows g) $ res
        , message_length  = getNRows g
        , codeword_length =  getNRows g + getNCols g
        }

mkLDPC_Code :: (MatrixLoader g, MatrixLoader h)
            => String
            -> (g -> U.Vector Bool -> IO (U.Vector Bool))
            -> (h -> Int -> U.Vector Double -> IO (U.Vector Bool))
            -> Code
mkLDPC_Code name encoder decoder = punctureTailOfCode $ Code ["ldpc/" ++ name ++ "/<matrix-name>/<max-rounds>[/.<truncation-size>]"]
     $ \ xs -> case xs of
                ["ldpc",nm,m,n] | nm == name && all isDigit n
                   -> fmap (: []) $ mkLDPC name m (read n) encoder decoder
                _  -> return []


-- Our version of atanh
{-# INLINE atanh' #-}
atanh' :: RealFloat a => a -> a
atanh' x | isInfinite y = signum x * 18.714973875118524
         | otherwise    = y
 where y = atanh x


newtype MinSum a = MinSum a
        deriving Show
{-# INLINE metric #-}
metric a b = min (abs a) (abs b) * signum a * signum b

instance (Ord a, Num a) => Semigroup (MinSum a) where
        {-# INLINE (<>) #-}
        (MinSum a) <> (MinSum b) = MinSum $ metric a b

data MinSum2 a = MinSum2 (MinSum a) Int (Maybe (MinSum a))
        deriving Show

instance  (Ord a, Num a) => Semigroup (MinSum2 a) where
        (MinSum2 a1@(MinSum v1) i1 b1) <> (MinSum2 a2@(MinSum v2) i2 b2)
                | abs v1 < abs v2 = MinSum2 (a1 <> a2) i1 (b1 <> Just a2)
                | otherwise       = MinSum2 (a1 <> a2) i2 (Just a1 <> b2 )

omit :: (Num a) => MinSum2 a -> Int -> a
omit (MinSum2 (MinSum a) i (Just (MinSum b))) j
        | i == j    = a         -- still needs sign adjust
        | otherwise = abs b     -- also needs sign adjust





