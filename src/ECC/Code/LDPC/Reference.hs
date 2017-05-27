{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Reference where

-- Reference implementation of LDPC

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Debug.Trace

type M a = Matrix a
type V a = U.Vector a

code :: Code
code = mkLDPC_Code "reference" encoder decoder

---------------------------------------------------------------------

encoder :: M Bool -> U.Vector Bool -> IO (U.Vector Bool)
encoder g v = return $ U.convert (getRow 1 (multStd' (rowVector $ U.convert v) (fmap fromBit (identity (nrows g)) <|> g)))

---------------------------------------------------------------------

decoder :: M Bool -> Int -> U.Vector Double -> IO (U.Vector Bool)
decoder a maxIterations orig_lam = return $ U.convert (ldpc a maxIterations (U.convert orig_lam))

toBit :: Bool -> Int
toBit False = 0
toBit True  = 1

fromBit :: Int -> Bool
fromBit 0 = False
fromBit 1 = True

multStd' :: M Bool -> M Bool -> M Bool
multStd' a b = fmap fromBit (multStd (fmap toBit a) (fmap toBit b))

ldpc :: M Bool -> Int -> V Double -> V Bool
ldpc a maxIterations orig_lam = U.map hard $ loop 0 orig_ne orig_lam
  where
    orig_ne :: M Double
    orig_ne = fmap (const 0) a

    loop :: Int -> M Double -> V Double -> V Double
    loop !n ne lam
        | U.all (== False) ans         = lam
        | n >= maxIterations           = orig_lam
        | otherwise                    = loop (n+1) ne' lam'
      where
        c_hat :: V Bool
        c_hat = U.map hard lam

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bool
        ans = U.convert $ getCol 1 (a `multStd'` colVector (U.convert c_hat))

        -- was bug here: V's start at index 0, not 1
        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a ! (m,n)
                then
                   -- was bug here: we need to cap atanh's answer
                    -2 * atanh' (product
                        [ tanh (- ((lam U.! (j-1) - ne ! (m,j)) / 2))
                        | j <- [1 .. U.length orig_lam]
                        , j /= n
                        , a ! (m,j)
                        ])
                else 0

        -- Was bug here: needed to add the orig_lam
        lam' :: V Double
        lam' = U.fromList [ U.foldr (+) (orig_lam U.! (j - 1)) (U.convert (getCol j ne'))
                          | j <- [1 .. U.length lam]
                          ]

