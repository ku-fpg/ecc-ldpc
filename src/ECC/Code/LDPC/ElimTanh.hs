{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.ElimTanh where

-- Implementation of LDPC, that approximates
-- 2 *tanh−1 ( PRODUCT [ tanh ( |Z(x[n′])| / 2)] )
--       with:
--  min |Z(x[n′])|
--       where n' eltOf N(m)\n (same for each equation)
--
-- See the paper, "Reduced-Complexity Decoding of LDPC Codes"
--  by J. Chen, A. Dholakia, E. Eleftheriou, M. Fossorier, and X.–Y. Hu


import Data.Bit
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import qualified Data.Vector as V

type M a = Matrix a
type V a = V.Vector a

code :: Code
code = Code ["ldpc/elimTanh/<matrix-name>/<max-rounds>[/<truncation-size>]"]
     $ \ xs -> return []
{-
Code ["ldpc/elimTanh/<matrix-name>/<max-rounds>[/<truncation-size>]"]
     $ \ xs -> case xs of
                ["ldpc","elimTanh",m,n]
                        | all isDigit n -> fmap (: []) $ mkLDPC m (read n)
                ["ldpc","elimTanh",m,n,t]
                        | all isDigit n
                       && all isDigit t -> fmap (: []) $ fmap (punctureTail (read t)) $ mkLDPC m (read n)
                _                       -> return []


mkLDPC :: String -> Int -> IO ECC
mkLDPC codeName maxI = do
   g :: G <- readAlist ("codes/" ++ codeName ++ ".G")   -- with G, we prepend the identity
   let full_g = identity (nrows g) <|> g
   h :: H <- readAlist ("codes/" ++ codeName ++ ".H")
   return $ ECC
        { name     = "ldpc/elimTanh/" ++ codeName ++ "/" ++ show maxI
        , encode   = return . V.toList . encoder full_g . V.fromList
        , decode   = return . (,True) . take (nrows full_g) . V.toList . ldpc maxI h . V.fromList
        , message_length  = nrows full_g
        , codeword_length =  ncols full_g
        }

type G = M Bit
type H = M Bit

encoder :: G -> V Bit -> V Bit
encoder g v = getRow 1 (multStd (rowVector v) g)

---------------------------------------------------------------------

ldpc :: Int -> M Bit -> V Double -> V Bit
ldpc maxIterations a orig_lam = fmap hard $ loop 0 orig_ne orig_lam
  where
    orig_ne :: M Double
    orig_ne = fmap (const 0) a

    loop :: Int -> M Double -> V Double -> V Double
    loop !n ne lam
        | V.all (== 0) ans             = lam
        | n >= maxIterations           = orig_lam
        | otherwise                    = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = fmap hard lam

        ans :: V Bit
        ans = getCol 1 (a `multStd` colVector c_hat)

        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a ! (m,n) == 1
                then
                    let bigVal = 100.0
                        (signProduct, absval) =
                            foldl ( \ (signProduct, curmin) x -> (signProduct /= ((signum x) == (-1)), min curmin (abs x)))
                                  (False, bigVal)
                                  [ lam V.! (j-1) - ne ! (m,j) | j <- [1 .. V.length orig_lam]
                                  , j /= n
                                  , a ! (m,j) == 1]
                    in if signProduct then negate absval else absval
                else 0

        lam' :: V Double
        lam' = V.fromList [ V.foldr (+) (orig_lam V.! (j - 1)) (getCol j ne')
                          | j <- [1 .. V.length lam]
                          ]
-}