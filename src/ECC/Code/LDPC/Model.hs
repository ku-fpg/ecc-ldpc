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
import qualified Data.Vector.Unboxed as U
import Data.Alist
import Debug.Trace
import Data.Semigroup
import qualified Data.BitMatrix.Word64 as BM64
import qualified Data.BitVector.Word64 as BV64
import qualified Data.Map as Map

import qualified ECC.Code.LDPC.Reference as Ref
import qualified Data.List.NonEmpty as NonEmpty

type M a = Matrix a
type V a = V.Vector a

code :: Code
code = mconcat
   [ mkLDPC_Code ("model-" ++ show n ++ "-" ++ show m) e d
   | (n,e) <- zip [0..] [ Ref.encoder, encoder1,encoder2,encoder3,encoder4,encoder1 `sameAs` encoder4]
   , (m,d) <- zip [0..] [ Ref.decoder
                        , decoder1
                        , decoder2
                        , decoder3
                        , decoder4
                        , decoder decoder5
                        , decoder decoder6
                        , decoder decoder7
                        , decoder decoder8
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

data Decoder m m_opt lam ne d ts i = Decoder
        { pre_a         :: m -> m_opt           -- ^ optimize the 'A' (/ H) matrix
        , pre_lambda    :: [Double] -> lam      -- ^ prepare the lambda
        , check_parity  :: m_opt -> lam -> Bool -- ^ check if the input meets parity
        , post_lambda   :: lam -> [Bit]         -- ^ output the result bits (with the parity bits truncated)
        , pre_ne        :: m_opt -> ne                  -- ^
        , comp_ne       :: forall ts . Semigroup ts => Share d ts i -> m_opt -> lam -> ne -> ne     -- ^
        , comp_lam      :: m_opt -> lam -> ne -> lam    -- ^ take the original lam and ne, and give back the new lam
        , share         :: Share d ts i
        }

data Share d ts i = Share
        { to_tanh       :: i -> d -> ts
        , from_tanh     :: ts -> i -> d
        }

sconcat' :: (Semigroup a) => [a] -> a
sconcat' = sconcat . NonEmpty.fromList

-- reference, recoded
decoder1 = decoder $ Decoder
        { pre_a        = id
        , pre_lambda   = V.fromList
        , check_parity = \ m_opt lam -> V.all (== 0) (getCol 1 (m_opt `multStd` colVector (fmap hard lam)))
        , post_lambda  = map hard . V.toList
        , pre_ne       = fmap (const 0)
        , comp_ne      = \ share m_opt lam ne -> matrix (nrows m_opt) (ncols m_opt) $ \ (m,n) ->
                if m_opt ! (m,n) == 1
                then
                    from_tanh share (sconcat'
                        [ to_tanh share j (lam V.! (j-1) - ne ! (m,j))
                        | j <- [1 .. ncols m_opt]
                        , m_opt ! (m,j) == 1
                        ]) n
                else 0
        , comp_lam     = \ m_opt orig_lam ne' ->
                V.fromList [ V.foldr (+) (orig_lam V.! (j - 1)) (getCol j ne')
                           | j <- [1 .. V.length orig_lam]
                           ]
        , share = share_tanh
        }

share_tanh :: (RealFloat d, Floating d, Fractional d, Eq i) => Share d [(i,d)] i
share_tanh = Share
                { to_tanh = \ j v -> [ (j,tanh (-v/2))]
                , from_tanh = \ jts n -> (-2) * atanh' (product [ t | (j,t) <- jts, j /= n ])
                }

data Tanh = Tanh Double

instance Monoid Tanh where
        mempty = Tanh 1
        mappend (Tanh a) (Tanh b) = Tanh $ a * b

-- using a finite map for the sparse array
decoder2 = decoder $ Decoder
        { pre_a        =  \ h ->( h
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                , V.fromList [ [ n | n <- [1..ncols h], h ! (m,n) == 1 ]
                                             | m <- [1..nrows h]
                                             ]
                                )
        , pre_lambda   = V.fromList
        , check_parity =  \ (m_opt,m,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (V.toList lam)))
        , post_lambda  =  map hard . V.toList
        , pre_ne       = \ (m_opt,_,_) -> Map.fromList
                                         [ ((m,n),0) | n <- [1..ncols m_opt], m <- [1..nrows m_opt], m_opt ! (m,n) == 1 ]
        , comp_ne      = \  share (m_opt,_,neighbors) lam ne -> Map.mapWithKey (\ (m,n) _ ->
                    from_tanh share (sconcat'
                        [ to_tanh share j (lam V.! (j-1) - ne Map.! (m,j))
                        | j <- neighbors V.! (m-1)
                        ]) n) ne
        , comp_lam     = \ (m_opt,_,_) orig_lam ne' ->
                V.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- Map.assocs ne' ]
        , share = share_tanh
        }

-- decoder3; pull out the call to to_tanh
decoder3 = decoder $ Decoder
        { pre_a        =  \ h -> ( h
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                , V.fromList [ [ n | n <- [1..ncols h], h ! (m,n) == 1 ]
                                             | m <- [1..nrows h]
                                             ]
                                )
        , pre_lambda   = V.fromList
        , check_parity =  \ (m_opt,m,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (V.toList lam)))
        , post_lambda  =  map hard . V.toList
        , pre_ne       = \ (m_opt,_,_) ->
                                Map.fromList [ ((m,n),0) | n <- [1..ncols m_opt], m <- [1..nrows m_opt], m_opt ! (m,n) == 1 ]
        , comp_ne      = \  share (m_opt,_,neighbors) lam ne ->
                let tmp = Map.mapWithKey (\ (_,n) v -> to_tanh share n (lam V.! (n-1) - v)) ne
                in Map.mapWithKey (\ (m,n) _ ->
                    from_tanh share (sconcat'
                        [ tmp Map.! (m,j)
                        | j <- neighbors V.! (m-1)
--                        , j /= n
                        ]) n) ne
        , comp_lam     = \ (m_opt,_,_) orig_lam ne' ->
                V.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- Map.assocs ne' ]
        , share = share_tanh
        }

decoder4 = decoder $ Decoder
        { pre_a        =  \ h -> ( h
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                , V.fromList [ [ n | n <- [1..ncols h], h ! (m,n) == 1 ]
                                             | m <- [1..nrows h]
                                             ]
                                )
        , pre_lambda   = V.fromList
        , check_parity =  \ (m_opt,m,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (V.toList lam)))
        , post_lambda  =  map hard . V.toList
        , pre_ne       = \ (m_opt,_,_) ->
                                Map.fromList [ ((m,n),0) | n <- [1..ncols m_opt], m <- [1..nrows m_opt], m_opt ! (m,n) == 1 ]
        , comp_ne      = \  share (m_opt,_,neighbors) lam ne ->
                let tanh_arr = V.fromList [ sconcat'
                                                [  to_tanh share j (lam V.! (j-1) - ne Map.! (m,j))
                                                | j <- neighbors V.! (m-1)
                                                ]
                                          | m <- [1..nrows m_opt]
                                          ] in
                Map.mapWithKey (\ (m,n) _ -> from_tanh share (tanh_arr V.! (m - 1)) n) ne
        , comp_lam     = \ (m_opt,_,_) orig_lam ne' ->
                V.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- Map.assocs ne' ]
        , share = share_tanh
        }


decoder5 = Decoder
        { pre_a        =  \ h ->
                                let vs = [ (m,n) | n <- [1..ncols h], m <- [1..nrows h], h ! (m,n) == 1 ] in
                                ( h
                                        -- The bit vector for the parity check
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                        -- all the left/right neibours
                                , V.fromList [ [ (i,j) | (i,(m',j))  <- [0..] `zip` vs, m == m' ]
                                             | m <- [1..nrows h]
                                             ]
                                        --
                                , V.fromList vs
                                )
        , pre_lambda   = V.fromList
        , check_parity =  \ (m_opt,m,_,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (V.toList lam)))
        , post_lambda  =  map hard . V.toList
        , pre_ne       = \ (m_opt,_,_,mns) -> V.map (const 0) mns
        , comp_ne      = \  share (m_opt,_,neighbors,mns) lam ne ->
                let tanh_arr = V.fromList [ sconcat'
                                                [ to_tanh share j ((lam V.! (j - 1)) - (ne V.! i))
                                                | (i,j) <- neighbors V.! (m-1)
                                                ]
                                          | m <- [1..nrows m_opt]
                                          ] in
                V.map (\ (m,n) -> from_tanh share (tanh_arr V.! (m - 1)) n) mns
        , comp_lam     = \ (m_opt,_,_,mns) orig_lam ne' ->
                V.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- V.toList mns `zip` V.toList ne' ]
        , share = share_tanh
        }

decoder6 = decoder5 { share = share_minsum }

-- Adding Nick's minsum2 optimization
decoder7 = Decoder
        { pre_a        =  \ h ->
                                let vs = [ (m,n) | n <- [1..ncols h], m <- [1..nrows h], h ! (m,n) == 1 ] in
                                ( h
                                        -- The bit vector for the parity check
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                        -- all the left/right neibours
                                , V.fromList [ [ (i,j) | (i,(m',j))  <- [0..] `zip` vs, m == m' ]
                                             | m <- [1..nrows h]
                                             ]
                                        --
                                , V.fromList vs
                                )
        , pre_lambda   = V.fromList
        , check_parity =  \ (m_opt,m,_,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (V.toList lam)))
        , post_lambda  =  map hard . V.toList
        , pre_ne       = \ (m_opt,_,_,mns) -> V.map (const 0) mns
        , comp_ne      = \  share (m_opt,_,neighbors,mns) lam ne ->
                -- the old way
                let tanh_arr = V.fromList [ sconcat'
                                                [ to_tanh share j ((lam V.! (j - 1)) - (ne V.! i))
                                                | (i,j) <- neighbors V.! (m-1)
                                                ]
                                          | m <- [1..nrows m_opt]
                                          ] in
                let ans1 = V.map (\ (m,n) -> from_tanh share (tanh_arr V.! (m - 1)) n) mns in

                -- The new way
                let interm_arr = V.zipWith (\ (_,n) v -> - ((lam V.! (n-1)) - v)) mns ne in

                let sign = V.accumulate (\ a b -> if b < 0 then not a else a)
                                   (V.generate (nrows m_opt) (const False))
                                   (V.zip (V.map (pred . fst) mns) interm_arr) in

                let val = V.accumulate (\ (b,c) a -> if abs a <= b
                                                     then (abs a,b)
                                                     else (b,min (abs a) c))
                                   (V.generate (nrows m_opt) (const (1/0,1/0)))     -- IEEE magic
                                   (V.zip (V.map (pred . fst) mns) interm_arr) in
                let ans2 = V.zipWith (\ (m,n) v ->
                                        let sgn = if sign V.! (m-1) == (v < 0) then 1 else -1 in
                                        let (a,b) = val V.! (m-1) in
                                        if a == abs v
                                        then (-0.75) * sgn * b
                                        else (-0.75) * sgn * a
                                     ) mns interm_arr in
                let signOf n | n < 0 = -1
                             | isNegativeZero n = -1
                             | otherwise = 1 in
                let comp = V.take 2
                         $ V.filter (\(a1,a2,_,b) -> not b)
                         $ (V.zipWith (\ a1 a2 -> (a1,a2, abs (a1 - a2), a1 == 0 || a1 == a2)
                                            ) ans1 ans2) in
                ans2
--                (traceShow comp ans1)
        , comp_lam     = \ (m_opt,_,_,mns) orig_lam ne' ->
                V.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- V.toList mns `zip` V.toList ne' ]
        , share = share_minsum
        }

-- Adding Unboxed arrays
decoder8 = Decoder
        { pre_a        =  \ h ->
                                let vs = [ (m,n) | n <- [1..ncols h], m <- [1..nrows h], h ! (m,n) == 1 ] in
                                ( h
                                        -- The bit vector for the parity check
                                , BM64.fromLists [[ h ! (m,n) | n <- [1..ncols h]] | m <- [1..nrows h]]
                                        -- all the left/right neibours
                                , V.fromList [ [ (i,j) | (i,(m',j))  <- [0..] `zip` vs, m == m' ]
                                             | m <- [1..nrows h]
                                             ]
                                        --
                                , U.fromList vs
                                )
        , pre_lambda   = U.fromList
        , check_parity =  \ (m_opt,m,_,_) lam -> not $ or $ BM64.parityMatVecMul m (BV64.fromList (fmap hard (U.toList lam)))
        , post_lambda  =  map hard . U.toList
        , pre_ne       = \ (m_opt,_,_,mns) -> U.map (const 0) mns
        , comp_ne      = \  share (m_opt,_,neighbors,mns) lam ne ->

                -- The new way
                -- Flat array of values
                let interm_arr = U.zipWith (\ (_,n) v -> - ((lam U.! (n-1)) - v)) mns ne in

                let sign = U.accumulate (\ a b -> if b < 0 then not a else a)
                                   (U.generate (nrows m_opt) (const False))
                                   (U.zip (U.map (pred . fst) mns) interm_arr) in

                let val = U.accumulate (\ (b,c) a -> if abs a <= b
                                                     then (abs a,b)
                                                     else (b,min (abs a) c))
                                   (U.generate (nrows m_opt) (const (1/0,1/0)))     -- IEEE magic
                                   (U.zip (U.map (pred . fst) mns) interm_arr) in
                let ans2 = U.zipWith (\ (m,_) v ->
                                        let sgn = if sign U.! (m-1) == (v < 0) then 1 else -1 in
                                        let (a,b) = val U.! (m-1) in
                                        if a == abs v
                                        then (-0.75) * sgn * b
                                        else (-0.75) * sgn * a
                                     ) mns interm_arr in
                ans2
--                (traceShow comp ans1)
        , comp_lam     = \ (m_opt,_,_,mns) orig_lam ne' ->
                U.accum (+) orig_lam [ (n-1,v) | ((_,n),v) <- U.toList mns `zip` U.toList ne' ]
        , share = share_minsum :: Share Double [(Int,Double)] Int -- ignored
        }



share_minsum :: (RealFloat d, Floating d, Fractional d, Eq i) => Share d [(i,d)] i
share_minsum = Share
                { to_tanh = \ j v -> [ (j,-v) ]
                , from_tanh = \ jts n -> (-0.75) * foldr1 (\ a b -> min (abs a) (abs b) * signum a * signum b)
                                                          [ t | (j,t) <- jts, j /= n ]

                }
{-

        = UnitMin2
        | Min2_1 Double
        | Min2_2 Double Double Double

instance Monoid Min2 where
        mempty                  = UnitMin2
        mappend a UnitMin       = a
        mappend UnitMin b       = b
        mappend (Min a) (Min b) = Min $ min (abs a) (abs b) * signum a * signum b
-}

getCol' j ne = map snd $ Map.toList $ Map.filterWithKey (\ (m,n) _ -> n == j) ne

-- V.toList (getRow 1 (multStd (rowVector $ V.fromList v) (identity (nrows g) <|> g)))

{-
ldpc :: M Bit -> Int -> V Double -> V Bit
ldpc a maxIterations orig_lam = fmap hard $ loop 0 orig_ne orig_lam
  where
    orig_ne :: M Double
    orig_ne = fmap (const 0) a
-}


decoder :: (Semigroup ts) => Decoder m m_opt lam ne d ts i -> m -> Int -> [Double] -> IO [Bit]
decoder dec a = \ !maxIterations inp -> do
      let orig_lam = pre_lambda dec inp

          loop !n ne lam
            | check_parity dec a_opt lam   = lam
            | n >= maxIterations           = orig_lam
            | otherwise                    =
                 let ne'  = {-# SCC ne' #-} comp_ne dec (share dec) a_opt lam      ne
                     lam' = {-# SCC lam' #-} comp_lam dec a_opt orig_lam ne'
                 in loop (n+1) ne' lam'

      return $ post_lambda dec $ loop 0 orig_ne orig_lam
   where
      a_opt   = pre_a dec a
      orig_ne = pre_ne dec a_opt

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
instance Same Bool where sameAs = defaultIsSame
instance Same Bit where sameAs = defaultIsSame

instance Same a => Same (IO a) where
    sameAs f1 f2 = liftM2 sameAs f1 f2

instance Same a => Same [a] where
    sameAs xs1 xs2 | length xs1 == length xs2 = zipWith sameAs xs1 xs2
                   | otherwise                = error ("lengths of lists are different" ++ show (length xs1,length xs2))




