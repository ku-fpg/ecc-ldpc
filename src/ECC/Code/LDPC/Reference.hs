{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Reference where

-- Reference implementation of LDPC

import Data.Bit
import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import qualified Data.Vector as V
import Data.Alist
import Debug.Trace

type M a = Matrix a
type V a = V.Vector a

code :: Code
code = mkLDPC_Code "reference" encoder decoder

---------------------------------------------------------------------

encoder :: M Bit -> [Bit] -> IO [Bit]
encoder g v = return $ V.toList (getRow 1 (multStd (rowVector $ V.fromList v) (identity (nrows g) <|> g)))

---------------------------------------------------------------------

decoder :: M Bit -> Int -> [Double] -> IO [Bit]
decoder a maxIterations orig_lam = return $ V.toList (ldpc a maxIterations (V.fromList orig_lam))

ldpc :: M Bit -> Int -> V Double -> V Bit
ldpc a maxIterations orig_lam = fmap hard $ loop 0 orig_ne orig_lam
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

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bit
        ans = getCol 1 (a `multStd` colVector c_hat)

        -- was bug here: V's start at index 0, not 1
        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a ! (m,n) == 1
                then
                   -- was bug here: we need to cap atanh's answer
                    -2 * atanh' (product
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

