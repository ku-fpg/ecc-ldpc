{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Reference where

-- Reference implementation of LDPC

import Data.Bit
import ECC.Types
import Data.Char (isDigit)

{-
import Data.Array.Matrix
import Data.Bit
import Data.Array
import ECC
-}

code :: Code
code = Code ["ldpc/reference/<matrix-name>/<max-rounds>"]
     $ \ xs -> case xs of
                        ["ldpc","reference",m,n]    -> fmap (: []) $ mkLDPC
                        _                          -> return []

mkLDPC :: IO ECC
mkLDPC = return ECC
        { name     = "ldpc/reference/"
        , encode   = return
        , decode   = return . (,True) . fmap mkBit . fmap (>= 0)
        , message_length  = 1
        , codeword_length = 1
        }


{-
import Data.Array.Matrix
import Data.Bit

encoder :: M Bit -> V Bit -> V Bit
encoder g v = getRowM $ rowM v `mm` g



decoder :: Int -> M Bit -> V Double -> V Bit
decoder = ldpc

ldpc :: forall d. (Floating d, Ord d) => Int -> M Bit -> V d -> V Bit
ldpc maxIterations a orig_lam = fmap hard $ loop 0 orig_ne orig_lam
  where
    orig_ne :: M d
    orig_ne = fmap (const 0) a

    loop :: Int -> M d -> V d -> V d
    loop !n ne lam
        | all (== 0) (elems ans)       = lam
        | n >= maxIterations           = orig_lam
        | otherwise                    = loop (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = fmap hard lam

        ans :: M Bit
        ans = a `mm` columnM c_hat

        ne' :: M d
        ne' = ne // [ ((m,n), -2 * atanh (product
                         [ tanh (- ((lam ! j - ne ! (m,j)) / 2))
                         | j <- indices lam
                         , j /= n
                         , a ! (m,j) == 1
                         ]))
                    | (m,n) <- indices ne
                    , a ! (m,n) == 1
                    ]

        lam' :: V d
        lam' = accum (+) orig_lam [ (n,a) | ((_,n),a) <- assocs ne' ]
-}