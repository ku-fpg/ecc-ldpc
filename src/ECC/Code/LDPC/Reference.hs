{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Reference where

-- Reference implementation of LDPC

import Data.Bit
import ECC.Types
import Data.Char (isDigit)
import Data.Matrix
import qualified Data.Vector as V
import Data.Alist

import Debug.Trace

type M a = Matrix a
type V a = V.Vector a

code :: Code
code = Code ["ldpc/reference/<matrix-name>/<max-rounds>"]
     $ \ xs -> case xs of
                        ["ldpc","reference",m,n]    -> fmap (: []) $ mkLDPC "moon.7.13" 200
                        _                          -> return []

mkLDPC :: String -> Int -> IO ECC
mkLDPC codeName maxI = do
   g :: G <- readAlist ("codes/" ++ codeName ++ ".G")
   h :: H <- readAlist ("codes/" ++ codeName ++ ".H")
   return $ ECC
        { name     = "ldpc/reference/"
        , encode   = return . V.toList . encoder g . V.fromList
        , decode   = return . (,True) . take (nrows g) . V.toList . ldpc 128 h . V.fromList
        , message_length  = nrows g
        , codeword_length =  ncols g
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
--        | traceShow ("loop",n,fmap round ne,lam) False   = undefined
        | V.all (== 0) ans             = lam
        | n >= maxIterations           = orig_lam
        | otherwise                    = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = fmap hard lam

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bit
        ans = getCol 1 (a `multStd` colVector c_hat)

        -- was bug here: V's start at index 0, not 1
        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a ! (m,n) == 1
                then
                    -2 * atanh (product
                        [ tanh (- ((lam V.! (j-1) - ne ! (m,j)) / 2))
                        | j <- [1 .. V.length orig_lam]
                        , j /= n
                        , a ! (m,j) == 1
                        ])
                else 0

        -- Was bug here: needed to add the orig_lam
        lam' :: V Double
        lam' = V.fromList [ V.foldr (+) (orig_lam V.! (j - 1)) (getCol j ne')
                          | j <- [1 .. V.length lam]
                          ]


crash = do
        h <- readAlist ("codes/moon.7.13.H")
        return $ ldpc 128 h v0

v0 = V.fromList $ [10.31704008699972,-5.762379657732054,10.71923587403389,8.69474407787769,12.403269631011147,-4.8140434138889425,1.096112697392211,-7.578778502331838,5.607400012378008,-7.687307594577467e-2,12.161984232264581,10.392572265154772,-5.1354846290270215,6.046205353391934,10.257982632106899,15.459395584777264,-5.998393431786397,12.594442761269947,1.387024485274938,5.832772771282147]
