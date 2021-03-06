{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections, KindSignatures, PolyKinds, ExistentialQuantification #-}
module ECC.Code.LDPC.Utils where

import ECC.Types
import ECC.Puncture
import Data.BitMatrix.Loader
import Data.Char (isDigit)
import Data.Monoid hiding ((<>))
import Data.Semigroup
import qualified Data.Vector.Unboxed as U
import Data.Ratio

import Control.Concurrent
import Data.IORef
import Control.Monad
import Data.List (foldl1')



-- | `mkLDPC` makes and LDPC EEC. The contact here is is that  
--   encode generates *at least* the (message_length - codeword_length),
--   decode can read upto the unpunctuated size, but generates
--   *at least* the message_length. 

mkLDPC :: (MatrixLoader g, MatrixLoader h)
       => String 
       -> String 
       -> Int
       -> Int
       -> Maybe (Ratio Int)
       -> (g -> Rate -> U.Vector Bool -> U.Vector Bool)
       -> (h -> IO (Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool))))
       -> IO (ECC IO)
mkLDPC prefix codeName maxI maxThreadCount optRate encoder decoder0 = do
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
   decoder1s <- replicateM maxThreadCount $ decoder0 h
   let decoder's = map (\f -> f rate maxI) decoder1s
   let unpuncture xs = U.take c_length xs `mappend` U.replicate (getNCols h - c_length) 0

   decoderInUse <- newIORef $ replicate maxThreadCount False

   foldl1' seq decoder's `seq` return $! ECC
        { name     = "ldpc/" ++ prefix ++ "/" ++ codeName ++ "/" ++ show maxI ++ "/" ++ show (numerator rate) ++ "/" ++  show (denominator rate)
        , encode   = \ inp -> pure (inp `mappend` U.take (c_length - m_length) (encoder' inp))
        , decode   = do
            tid <- myThreadId
            let tidN = read $ drop (length "ThreadId" + 1) (show tid) :: Int
            -- print =<< fmap (\x -> ("myThreadId", x)) myThreadId
            -- print =<< threadCapability =<< myThreadId

            pure $ \ inp -> do
              r0 <- (decoder's !! (tidN `rem` maxThreadCount)) (unpuncture inp)
              pure $ case r0 of
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
mkLDPC_Code name encoder decoder = Code ["ldpc/" ++ name ++ "/<matrix-name>/<max-rounds>[/codeword/message]"] (pure ()) (const (pure ()))
     $ \ vars xs -> case xs of
                ["ldpc",nm,m,n,x,y] | nm == name && all isDigit n && all isDigit x && all isDigit y
                   -> decoder `seq` fmap (: []) $ mkLDPC name m (read n) 1 (Just (read x % read y)) encoder (pure . decoder')
                ["ldpc",nm,m,n] | nm == name && all isDigit n
                   -> decoder `seq` fmap (: []) $ mkLDPC name m (read n) 1 Nothing encoder (pure . decoder')
                _  -> return []
      where decoder' h = \x y z -> pure $ decoder h x y z

mkLDPC_CodeIO :: (MatrixLoader g, MatrixLoader h)
            => String
            -> Int
            -> (g -> Rate -> U.Vector Bool -> U.Vector Bool)
            -> (vars -> h -> IO (Rate -> Int -> U.Vector Double -> IO (Maybe (U.Vector Bool))))
            -> IO vars
            -> (vars -> IO ())
            -> Code
mkLDPC_CodeIO name maxThreadCount encoder decoder initialize finalize =
  Code ["ldpc/" ++ name ++ "/<matrix-name>/<max-rounds>[/codeword/message]"]
       initialize
       finalize
       $ \ vars xs -> case xs of
                        ["ldpc",nm,m,n,x,y] | nm == name && all isDigit n && all isDigit x && all isDigit y
                           -> fmap (: []) $ mkLDPC name m (read n) maxThreadCount (Just (read x % read y)) encoder (decoder vars)
                        ["ldpc",nm,m,n] | nm == name && all isDigit n 
                           -> fmap (: []) $ mkLDPC name m (read n) maxThreadCount Nothing encoder (decoder vars)
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

{-# INLINE metricDouble #-}
metricDouble :: Double -> Double -> Double
metricDouble a b = min (abs a) (abs b) * signum a * signum b

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



