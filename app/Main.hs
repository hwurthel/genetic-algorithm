module Main where

import StopCondition
import Evolution
import Protein
import IO

import Data.List (sortOn)

main :: IO ()
main = do
    let stepOfEvolution = \x xs -> x >>= crossover >>= mutation >>= flip computeLambda xs >>= selection
        process n sc (pop_x, all_x, best_x) = do
            -- Делаем шаг эволюции
            pop_x'  <- stepOfEvolution pop_x all_x

            -- Печатаем новую популяцию в файл
            writeInFile ("Step " <> show n <> "\n") pop_x'
            
            -- Дополняем множество всех особей новой популяцией 
            all_x' <- return $ all_x <> pop_x'

            -- Ищем лучшую особь в новой популяции и сравниваем
            -- её с лучшей среди особей во всех популяциях
            let (best_x', sc') = if maximum pop_x' > best_x 
                then (maximum pop_x', 0)
                else (best_x, sc + 1) 

            -- Проверяем условия остановки.
            -- В зависимости от результата возвращаем
            -- все особи или продолжаем работу
            if isStop sc' then return all_x'
            else process (n + 1) sc' (return pop_x', all_x', best_x')

    res <- process 1 0 (generatePopulation, [], tmpProtein)
    writeInFile "--------------\n--------------\n" []
    writeInFile "Result:\n" (reverse $ sortOn lambda res)
    return ()