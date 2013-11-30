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

type M a = Matrix a
type V a = V.Vector a

{-
import Data.Array.Matrix
import Data.Bit
import Data.Array
import ECC
-}

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
        , encode   = return . encoder g
        , decode   = return . (,True) . fmap mkBit . fmap (>= 0)
        , message_length  = nrows g
        , codeword_length =  ncols g
        }

type G = M Bit
type H = M Bit

encoder :: G -> [Bit] -> [Bit]
encoder g v = V.toList (getRow 1 (multStd (rowVector (V.fromList v)) g))

{-
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