{-# LANGUAGE PatternSynonyms #-}

module Main (main) where

import Data.Matrix
import System.Environment

import System.FilePath.Posix (splitExtension)

import Data.Matrix.Matlab

pattern Extension = ".qx"

main :: IO ()
main = do
  args <- getArgs

  case args of
    [fileName]
      | (_, Extension) <- splitExtension fileName -> do
          fileContents <- readFile fileName

          let Matlab mat = read fileContents :: Matlab Integer

          print (Matlab (fmap (^2) mat))
      | otherwise ->
          error ("Wrong extension (expected " ++ Extension ++ ")")

    _ -> error "Wrong number of arguments. Expected 1."

