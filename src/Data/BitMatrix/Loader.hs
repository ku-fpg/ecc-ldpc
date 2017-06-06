{-# LANGUAGE FlexibleInstances, GADTs #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE FlexibleContexts #-}
module Data.BitMatrix.Loader where

--import qualified Data.ByteString.Lazy as LBS
import System.Directory
import Data.BitMatrix.Alist
import Data.BitMatrix.Matlab
import Data.Matrix.QuasiCyclic
import Data.Matrix
import System.IO
import Data.Char(chr)
import qualified Codec.Compression.GZip as GZip
import Control.Monad

import Paths_ecc_ldpc

-- read a matrix from a file.
-- We have different types of matrix formats
-- * matrix (.m)
-- * alist (.alist)
-- * quazi-cyclic (.q<N>)
-- We also support
-- * compression (.gz)
--
-- So we support multi-unpacking, based on the name.

class MatrixLoader m where
        getMatrix ::  LoaderMatrix -> m
        getNRows   :: m -> Int                  -- number of rows in the *expanded* matrix
        getNCols   :: m -> Int                  -- number of cols in the *expanded* matrix

instance MatrixLoader (Matrix Bool) where
     getMatrix (LoaderAline m)  = m
     getMatrix (LoaderMatlab m) = m
     getMatrix (LoaderQC m)     = toBitMatrix m
     getNRows = nrows
     getNCols = ncols

instance MatrixLoader (QuasiCyclic Integer) where
     getMatrix (LoaderAline m)  = error "can not load aline as QuasiCyclic"
     getMatrix (LoaderMatlab m) = error "can not load matlab as QuasiCyclic"
     getMatrix (LoaderQC m)     = m
     getNRows (QuasiCyclic n a) = n * nrows a
     getNCols (QuasiCyclic n a) = n * ncols a


-- The closed internally supported Datatypes, with their efficent representations.
data LoaderMatrix where
        LoaderAline  :: Matrix Bool         -> LoaderMatrix
        LoaderMatlab :: Matrix Bool         -> LoaderMatrix
        LoaderQC     :: QuasiCyclic Integer -> LoaderMatrix
  deriving Show

-- deriving instance Show LoaderMatrix

loaders :: [(String,FilePath -> IO LoaderMatrix)]
loaders = [("alist",    fmap (LoaderAline . unAlist . read)     . readFile)
          ,("m",        fmap (LoaderMatlab . unMatlab . read)   . readFile)
          ,("q",        fmap (LoaderQC . read)                  . readFile)
          ]

loadMatrix :: MatrixLoader m => FilePath -> IO m
loadMatrix filePath = do
        dat <- getDataDir
        print dat
        let prefixes = ["codes",dat ++ "/codes"]
        filePaths <- sequence [ do ok <- doesFileExist file
                                   return (ok,loader file)
                              | (suffix,loader) <- loaders
                              , prefix <- prefixes
                              , let file = prefix ++ "/" ++ filePath ++ "." ++ suffix
                              ]
        case [ loader | (True,loader) <- filePaths ] of
                []         -> error $ "can not find any matrix files in " ++ show filePath
                (loader:_) -> loader >>= return . getMatrix

