-- | Glue (row) consecutive unsigned 64-bit integers into 128-bit integers.
{-# LANGUAGE PatternSynonyms #-}

module Main (main) where

import Data.Word (Word64)

import Data.Matrix
import System.Environment

import System.FilePath.Posix (splitExtension)

import Data.Matrix.Matlab

import Data.Bits

pattern Extension = ".q64"

rowGlue :: [Word64] -> [Integer]
rowGlue (x:y:rest) = ((fromIntegral x `shiftL` 64) .|. fromIntegral y) : rowGlue rest
rowGlue []         = []
rowGlue _          = error "Odd number of items in row."

main :: IO ()
main = do
  args <- getArgs

  case args of
    [fileName]
      | (_, Extension) <- splitExtension fileName -> do
          fileContents <- readFile fileName

          let Matlab mat = read fileContents :: Matlab Word64

              mat'Lists  = rowGlue <$> toLists mat

          print (Matlab (fromLists mat'Lists))

      | otherwise ->
          error ("Wrong extension (expected " ++ Extension ++ ")")

    _ -> error "Wrong number of arguments. Expected 1."


