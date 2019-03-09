module Main where

import StopCondition
import Evolution
import Protein
import IO

main :: IO ()
main = do
    let stepOfEvolution = \x -> x >>= crossover >>= mutation >>= computeLambda >>= selection
        process n (x, xs) = do
            x'  <- stepOfEvolution x
            xs' <- return $ xs <> x'
            writeInFile ("Step " <> show n <> "\n") x'
            if isStop n then return xs'
            else process (n + 1) (return x', xs') 
    res <- process 1 (generatePopulation, [])
    writeInFile "--------------\n--------------\n" []
    writeInFile "Result:\n" (sortPopulation res)
    return ()