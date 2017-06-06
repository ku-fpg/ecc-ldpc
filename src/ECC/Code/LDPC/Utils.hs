{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Utils where

import ECC.Types
import ECC.Puncture
import Data.BitMatrix.Loader
import Data.Char (isDigit)
import Data.Monoid hiding ((<>))
import Data.Semigroup
import qualified Data.Vector.Unboxed as U
import Data.Ratio



-- | `mkLDPC` makes and LDPC EEC. The contact here is is that  
--   encode generates *at least* the (message_length - codeword_length),
--   decode can read upto the unpunctuated size, but generates
--   *at least* the message_length. 

mkLDPC :: (MatrixLoader g, MatrixLoader h)
       => String 
       -> String 
       -> Int
       -> Maybe (Ratio Int)
       -> (g -> Rate -> U.Vector Bool -> U.Vector Bool)
       -> (h -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool))
       -> IO (ECC IO)
mkLDPC prefix codeName maxI optRate encoder decoder = do
   g <- loadMatrix (codeName ++ "/G") -- with G, we prepend the identity
   print (getNRows g, getNCols g)
   -- #rows are the message size, #cols + #rows are the (unpunctuated) message
   let m_length        = getNRows g
   h <- loadMatrix (codeName ++ "/H")
   -- #rows are number of paritiy checks, #cols are the (unpunctuated) message
   print (getNRows h, getNCols h)
   if getNRows g + getNCols g /= getNCols h then error ("bad code size match" ++ show (getNRows g + getNCols g, getNCols h))
                                            else return ()
   let unpunc_c_length = getNCols h
   let rate = case optRate of
               Just r -> r
               Nothing -> m_length % unpunc_c_length
   print $ rate
   let c_length = (m_length * denominator rate) `div` numerator rate
   print $ c_length
   let encoder' = encoder g rate
   let decoder' = decoder h rate maxI
   let unpuncture xs = U.take c_length xs `mappend` U.replicate (getNCols h - c_length) 0
   return $ ECC
        { name     = "ldpc/" ++ prefix ++ "/" ++ codeName ++ "/" ++ show maxI ++ "/" ++ show (numerator rate) ++ "/" ++  show (denominator rate)
        , encode   = \ inp -> pure (inp `mappend` U.take (c_length - m_length) (encoder' inp))
        , decode   = \ inp -> pure 
                            $ case decoder' (unpuncture inp) of
                                Nothing  -> (U.take m_length  $ U.map hard inp, False)
                                Just  r  -> (U.take m_length $ r, True)
        , message_length = m_length
        , codeword_length = c_length
        }

mkLDPC_Code :: (MatrixLoader g, MatrixLoader h)
            => String
            -> (g -> Rate -> U.Vector Bool -> U.Vector Bool)
            -> (h -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool))
            -> Code
mkLDPC_Code name encoder decoder = Code ["ldpc/" ++ name ++ "/<matrix-name>/<max-rounds>[/codeword/message]"]
     $ \ xs -> case xs of
                ["ldpc",nm,m,n,x,y] | nm == name && all isDigit n && all isDigit x && all isDigit y
                   -> fmap (: []) $ mkLDPC name m (read n) (Just (read x % read y)) encoder decoder
                ["ldpc",nm,m,n] | nm == name && all isDigit n 
                   -> fmap (: []) $ mkLDPC name m (read n) Nothing encoder decoder
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





