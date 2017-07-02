{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ScopedTypeVariables #-}

module ECC.Code.LDPC.GPU.Reference where

import Prelude hiding (map, Num, Eq, (>), (<), (>=), (<=), (==), (/=), zipWith, all, (||), (++), fst, snd, (!!), not, replicate, (&&), Fractional, Floating, RealFloat, filter, zip, any, even, sum, or, and, odd)
import qualified Prelude as P

import ECC.Code.LDPC.Utils
import ECC.Types
import ECC.Code.LDPC.Fast.Encoder (encoder)
import ECC.Puncture
import Data.Char (isDigit)
import qualified Data.Matrix as M
import Data.Bit
import qualified Data.Vector.Unboxed as U
import qualified Data.Vector as V
import Debug.Trace

import Data.Array.Accelerate as A
import Data.Array.Accelerate.Debug
import Data.Array.Accelerate.IO
import Data.Array.Accelerate.LLVM.PTX

import Control.Arrow (first)

type V a = Array DIM1 a
type M a = Array DIM2 a

code :: Code
code = compiledLdpc `seq` mkLDPC_Code "gpu-reference" encoder decoder

---------------------------------------------------------------------

decoder :: M.Matrix Bool -> Rate -> Int -> U.Vector Double -> Maybe (U.Vector Bool)
decoder a0 rate maxIterations orig_lam =
  Just $
  U.map word8ToBool $
  U.convert $
  toVectors $
  compiledLdpc
    (a
    ,run $ unit $ lift maxIterations
    ,(fromVectors (Z :. U.length orig_lam) (U.convert orig_lam))
    )
  where
    a = fromVectors (Z :. M.nrows a0 :. M.ncols a0) (U.convert (fmap boolToWord8 (M.getMatrixAsVector a0)))

word8ToBool :: Word8 -> Bool
word8ToBool 0 = False
word8ToBool 1 = True
word8ToBool n = error $ "word8ToBool: " P.++ show n

boolToWord8 :: Bool -> Word8
boolToWord8 False = 0
boolToWord8 True  = 1

hard' :: Exp Double -> Exp Bool
hard' = (> 0)

{-# INLINE atanh'' #-}
atanh'' :: (RealFloat a, Eq a, Fractional a, Floating a) => Exp a -> Exp a
atanh'' x =
  cond (x == 1 || x == -1 || atanhX >= inf || atanhX <= negInf)
       (signum x * 18.714973875118524)
       atanhX
  where
    atanhX = atanh x

inf, negInf :: (RealFloat a, Floating a) => Exp a
inf    =  1/0
negInf = -1/0

compiledLdpc :: (M Bool, Scalar Int, V Double) -> V Bool
compiledLdpc = run1 ldpc

ldpc :: Acc (M Bool, Scalar Int, V Double) -> Acc (V Bool)
-- ldpc a0 maxIterations orig_lam = map hard' $ loop (unit 0) orig_ne orig_lam
ldpc t0 = map hard' $ loop (unit 0) orig_ne orig_lam
  where
    (a0, maxIterations, orig_lam) = unlift t0 :: (Acc (M Bool), Acc (Scalar Int), Acc (V Double))

    orig_ne :: Acc (M Double)
    orig_ne = fill (shape a0) 0

    loop :: Acc (Scalar Int) -> Acc (M Double) -> Acc (V Double) -> Acc (V Double)
    loop n ne lam =
      case unlift $
             awhile loopCond
                    loopBody
                    (lift (n, ne, lam)) :: (Acc (Scalar Int), Acc (M Double), Acc (V Double))
           of
        (n', ne', lam') ->
          acond (the n' >= the maxIterations)
                orig_lam
                lam'

    loopCond :: Acc (Scalar Int, M Double, V Double) -> Acc (Scalar Bool)
    loopCond t = unit $ (the n < the maxIterations) && the continue
      where

        continue :: Acc (Scalar Bool)
        continue =
          or $
          imap
            (\ix _ ->
              let r = unindex1 $ unlift ix :: Exp Int
              in
              odd . snd $
              while
                ((< a0ColCount) . fst)
                (\p ->
                  let (c, v) = unlift p :: (Exp Int, Exp Int)
                  in
                  lift
                    (c+1
                    ,cond (a0 ! lift (Z :. r :. c) && c_hat ! lift (Z :. c))
                          (v+1)
                          v
                    )
                )
                (lift (0::Int,0::Int)))
            a0Rows

        -- -- Parity of the population count of each row
        -- ans :: Acc (V Bool)
        -- ans =
        --   generate
        --     (lift (Z :. a0RowCount))
        --     (\ix ->
        --       let r = unindex1 $ unlift ix :: Exp Int
        --       in
        --       even . snd $
        --       while
        --         ((< a0ColCount) . fst)
        --         (\p ->
        --           let (c, v) = unlift p :: (Exp Int, Exp Int)
        --           in
        --           lift
        --             (c+1
        --             ,cond (a0 ! lift (Z :. r :. c) && c_hat ! lift (Z :. c))
        --                   v
        --                   (v+1)
        --             )
        --         )
        --         (lift (0::Int,0::Int)))


        c_hat :: Acc (V Bool)
        c_hat = map hard' lam

        (n, _, lam) = unlift t :: (Acc (Scalar Int), Acc (M Double), Acc (V Double))

    loopBody :: Acc (Scalar Int, M Double, V Double) -> Acc (Scalar Int, M Double, V Double)
    loopBody t =
      lift (unit (the iters + 1), ne', lam')
      where
        (iters, ne, lam) = unlift t :: (Acc (Scalar Int), Acc (M Double), Acc (V Double))

        ne' :: Acc (M Double)
        ne' = generate (shape ne) (\ ix ->
          let (m, n) = unlift $ unindex2 ix :: (Exp Int, Exp Int)
          in
          cond (not (a0 ! ix))
               0
               (
               let wh =
                        (while
                     ((< colCount) . fst)
                     (\ p ->
                         let (j, v) = unlift p :: (Exp Int, Exp Double)
                         in
                         cond
                           (j == n || not (a0 ! (lift (Z :. m :. j))))
                           (lift (j+1, v))
                           (lift
                             (j+1
                             -- tanh (- ((lam U.! (j-1) - ne ! (m,j)) / 2))
                             ,v*tanh (- ((lam ! lift (Z :. j)) - (ne ! lift (Z :. m :. j)))/2)
                             )))
                     (lift (0 :: Exp Int, 1 :: Exp Double)
                        :: Exp (Int, Double)))
               in
               let (_, r) = unlift wh :: (Exp Int, Exp Double)
               in
               -2 * atanh'' r
          ))

        lam' :: Acc (V Double)
        lam' = zipWith (+) orig_lam (fold1 (+) (transpose ne'))

    (a0RowCount, a0ColCount) = unlift $ unindex2 (shape a0) :: (Exp Int, Exp Int)
    a0Rows = enumFromN (lift (Z :. a0RowCount)) 0 :: Acc (V Int)
    colCount = unindex1 (shape orig_lam) :: Exp Int


matrixVecMult :: Acc (M Bool) -> Acc (V Bool) -> Acc (V Bool)
matrixVecMult m v =
  fold1 (/=) $ zipWith (&&) m (replicate (lift (Z :. All :. colCount)) v)
  where
    (_, colCount) = unlift $ unindex2 (shape m) :: (Exp Int, Exp Int)
