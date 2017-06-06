{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Reference.Orig where

-- Reference implementation of LDPC

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Debug.Trace

type M a = Matrix a
type V a = U.Vector a

code :: Code
code = mkLDPC_Code "reference" encoder decoder

---------------------------------------------------------------------

encoder :: M Bool -> Rate -> U.Vector Bool -> U.Vector Bool
encoder g _ v = U.map toBool $ U.convert (getRow 1 (multStd (rowVector $ U.convert (U.map fromBool v)) ((fmap fromBool g))))

---------------------------------------------------------------------

decoder :: M Bool -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a rate maxIterations orig_lam = Just $ U.convert (ldpc a maxIterations (U.convert orig_lam))

toBit :: Num a => Bool -> a
toBit False = 0
toBit True  = 1

fromBit :: (Num a, Eq a, Show a) => a -> Bool
fromBit 0 = False
fromBit 1 = True
fromBit n = error $ "fromBit: " ++ show n

multStd' :: M Bool -> M Bool -> M Bool
multStd' a b = fmap fromBit $ multStd (fmap toBit a) (fmap toBit b)

toBits :: M Bool -> M Bit
toBits = fmap toBit

fromBits :: M Bit -> M Bool
fromBits = fmap fromBit

ldpc :: M Bool -> Int -> V Double -> V Bool
ldpc a0 maxIterations orig_lam = U.map hard $ loop 0 orig_ne orig_lam
  where
    a :: M Bit
    a = toBits a0

    orig_ne :: M Double
    orig_ne = fmap (const 0) a

    loop :: Int -> M Double -> V Double -> V Double
    loop !n ne lam
        | U.all (== False) ans     = lam
        | n >= maxIterations       = orig_lam
        | otherwise                = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = U.map (toBit . hard) lam

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bool
        ans = U.convert $ fmap fromBit $ getCol 1 (a `multStd` colVector (U.convert c_hat))

        -- was bug here: V's start at index 0, not 1
        ne' :: M Double
        ne' = matrix (nrows orig_ne) (ncols orig_ne) $ \ (m,n) ->
                if a ! (m,n) == 1
                then
                   -- was bug here: we need to cap atanh's answer
                    -2 * atanh' (product
                        [ tanh (- ((lam U.! (j-1) - ne ! (m,j)) / 2))
                        | j <- [1 .. U.length orig_lam]
                        , j /= n
                        , a ! (m,j) == 1
                        ])
                else 0

        -- Was bug here: needed to add the orig_lam
        lam' :: V Double
        lam' = U.fromList [ U.foldr (+) (orig_lam U.! (j - 1)) (U.convert (getCol j ne'))
                          | j <- [1 .. U.length lam]
                          ]

