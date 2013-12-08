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
   | (n,e) <- zip [1..] [encoder1,encoder2,encoder3,encoder4]
   , (m,d) <- zip [1..] [ \ (_ :: Matrix Bit) _ vs -> (sumBits (map hard vs) == 0) `seq` return (map hard vs)
                        , Ref.decoder
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

encoder4 = encoder (Encoder { pre_matrix     = \ g -> BM64.fromLists [[ g ! (m,n) | n <- [1..nrows g]] | m <- [1..ncols g]]
                            , prop_input     = \ vs -> (vs,BV64.fromList vs)
                            , mat_mul        = \ (vs,v) m -> (vs,BM64.vecMatMul v m)
                            , post_output    = \ (vs,v) -> vs ++ BV64.toList v
                            })



encoder :: Encoder m m_opt inp_opt out_opt -> m -> [Bit] -> IO [Bit]
encoder enc g = let m_opt = pre_matrix enc g
                in \ v -> return $ post_output enc (mat_mul enc (prop_input enc v) m_opt)

-- V.toList (getRow 1 (multStd (rowVector $ V.fromList v) (identity (nrows g) <|> g)))

decoder :: M Bit -> Int -> [Double] -> IO [Bit]
decoder a maxIterations orig_lam = return $ V.toList (ldpc a maxIterations (V.fromList orig_lam))

ldpc :: M Bit -> Int -> V Double -> V Bit
ldpc a maxIterations orig_lam = error "" -- fmap hard $ loop 0 orig_ne orig_lam $ error ""
