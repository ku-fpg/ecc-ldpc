{-# LANGUAGE BangPatterns, RankNTypes, GeneralizedNewtypeDeriving, ScopedTypeVariables #-}
{-# LANGUAGE TupleSections #-}
module ECC.Code.LDPC.Model where

-- Reference implementation of LDPC

import Data.Bit
import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Puncture
import Data.Char (isDigit)
import Data.Matrix
import Control.Monad
import qualified Data.Vector as V
import Data.Alist
import Debug.Trace
import Data.Monoid
import qualified Data.BitMatrix.Word64 as BM64
import qualified Data.BitVector.Word64 as BV64

import qualified ECC.Code.LDPC.Reference as Ref

type M a = Matrix a
type V a = V.Vector a

code :: Code
code = mconcat
   [ mkLDPC_Code ("model-" ++ show n ++ "-" ++ show m) e d
   | (n,e) <- zip [0..] [ Ref.encoder, encoder1,encoder2,encoder3,encoder4,encoder1 `sameAs` encoder4]
   , (m,d) <- zip [0..] [ Ref.decoder
                        , decoder1
                        , \ (_ :: Matrix Bit) _ vs -> (sumBits (map hard vs) == 0) `seq` return (map hard vs)
                        ]
   ]


data Encoder m m_opt inp_opt out_opt = Encoder
        { pre_matrix     :: m       -> m_opt                    -- prepare the matrix
        , prop_input     :: [Bit]   -> inp_opt                  -- process input
        , mat_mul        :: inp_opt -> m_opt -> out_opt         -- multiply
        , post_output    :: out_opt -> [Bit]                    --
        }

encoder1 = encoder (Encoder { pre_matrix     = \ g -> identity (nrows g) <|> g
                            , prop_input     = rowVector . V.fromList
                            , mat_mul        = multStd
                            , post_output    = V.toList . getRow 1
                            })

encoder2 = encoder (Encoder { pre_matrix     = \ g -> g
                            , prop_input     = \ vs -> (vs,rowVector (V.fromList vs))
                            , mat_mul        = \ (vs,v) m -> (vs,multStd v m)
                            , post_output    = \ (vs,v) -> vs ++ V.toList (getRow 1 v)
                            })

-- dumbest possible encoder, fails to include parity
encoder3 = encoder (Encoder { pre_matrix     = \ (m :: Matrix Bit) -> ncols m + nrows m
                            , prop_input     = \ vs -> vs
                            , mat_mul        = \ vs m -> take m (vs ++ repeat 0)
                            , post_output    = \ vs -> vs
                            })

encoder4 = encoder (Encoder { pre_matrix     = \ g -> BM64.fromLists [[ g ! (m,n) | n <- [1..ncols g]] | m <- [1..nrows g]]
                            , prop_input     = \ vs -> (vs,BV64.fromList vs)
                            , mat_mul        = \ (vs,v) m -> (vs,BM64.vecMatMul v m)
                            , post_output    = \ (vs,v) -> vs ++ BV64.toList v
                            })


encoder :: Encoder m m_opt inp_opt out_opt -> m -> [Bit] -> IO [Bit]
encoder enc g = let m_opt = pre_matrix enc g
                in \ v -> return $ post_output enc (mat_mul enc (prop_input enc v) m_opt)

----------------------------------------------------

data Decoder m m_opt lam ne = Decoder
        { pre_a         :: m -> m_opt           -- ^ optimize the 'A' (/ H) matrix
        , pre_lambda    :: [Double] -> lam      -- ^ prepare the lambda
        , check_parity  :: m_opt -> lam -> Bool -- ^ check if the input meets parity
        , post_lambda   :: lam -> [Bit]         -- ^ output the result bits (with the parity bits truncated)
        , pre_ne        :: m_opt -> ne                  -- ^
        , comp_ne       :: m_opt -> lam -> ne -> ne     -- ^
        , comp_lam      :: m_opt -> lam -> ne -> lam    -- ^ take the original lam and ne, and give back the new lam
        }


decoder1 = decoder $ Decoder
        { pre_a        = id
        , pre_lambda   = V.fromList
        , check_parity = \ m_opt lam -> V.all (== 0) (getCol 1 (m_opt `multStd` colVector (fmap hard lam)))
        , post_lambda  = map hard . V.toList
        , pre_ne       = fmap (const 0)
        , comp_ne      = \ m_opt lam ne -> matrix (nrows m_opt) (ncols m_opt) $ \ (m,n) ->
                if m_opt ! (m,n) == 1
                then
                    -2 * atanh' (product
                        [ tanh (- ((lam V.! (j-1) - ne ! (m,j)) / 2))
                        | j <- [1 .. ncols m_opt]
                        , j /= n
                        , m_opt ! (m,j) == 1
                        ])
                else 0
        , comp_lam     = \ m_opt orig_lam ne' ->
                V.fromList [ V.foldr (+) (orig_lam V.! (j - 1)) (getCol j ne')
                           | j <- [1 .. V.length orig_lam]
                           ]
        }

-- V.toList (getRow 1 (multStd (rowVector $ V.fromList v) (identity (nrows g) <|> g)))

{-
ldpc :: M Bit -> Int -> V Double -> V Bit
ldpc a maxIterations orig_lam = fmap hard $ loop 0 orig_ne orig_lam
  where
    orig_ne :: M Double
    orig_ne = fmap (const 0) a
-}

decoder :: Decoder m m_opt lam ne -> m -> Int -> [Double] -> IO [Bit]
decoder dec a = \ !maxIterations inp -> do
      let orig_lam = pre_lambda dec inp

          loop !n ne lam
            | check_parity dec a_opt lam   = lam
            | n >= maxIterations           = orig_lam
            | otherwise                    =
                 let ne'  = comp_ne  dec a_opt lam      ne
                     lam' = comp_lam dec a_opt orig_lam ne'
                 in loop (n+1) ne' lam'

      return $ post_lambda dec $ loop 0 (pre_ne dec a_opt) orig_lam
   where
      a_opt = pre_a dec a

{-
        | otherwise                    = loop maxIterations (n+1) ne' lam'
      where
        c_hat :: V Bit
        c_hat = fmap hard lam

        -- was bug here: needed to getCol, not getRow (which was a single element)
        ans :: V Bit
        ans =

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
-}


-- maxIterations orig_lam = return $ V.toList (ldpc a maxIterations (V.fromList orig_lam))

ldpc :: M Bit -> Int -> V Double -> V Bit
ldpc a maxIterations orig_lam = error "" -- fmap hard $ loop 0 orig_ne orig_lam $ error ""

class Same a where
  sameAs :: Same a => a -> a -> a

{-
sameAs x y = case isSame x y of
               Left msg -> error $ "sameAs failed: " ++ msg
               Right a -> a
-}

defaultIsSame :: Eq a => a -> a -> a
defaultIsSame a b | a == b    = a
                  | otherwise = error "equality failure"

instance (Eq a, Same b) => Same (a -> b) where
    sameAs f1 f2 x = f1 x `sameAs` f2 x

instance Same Int where sameAs = defaultIsSame
instance Same Bit where sameAs = defaultIsSame

instance Same a => Same (IO a) where
    sameAs f1 f2 = liftM2 sameAs f1 f2

instance Same a => Same [a] where
    sameAs xs1 xs2 | length xs1 == length xs2 = zipWith sameAs xs1 xs2
                   | otherwise                = error ("lengths of lists are different" ++ show (length xs1,length xs2))




