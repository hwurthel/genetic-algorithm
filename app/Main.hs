module Main where

import StopCondition
import Evolution
import Protein
import InputOutput
import InsertMolecule

import System.Environment
import Data.List
import Data.Maybe

insertMolecule :: IO ()
insertMolecule = do
    args <- getArgs
    let zmatr_name = args !! 0
        mol_name   = args !! 1
        n = read $ args !! 2
        s = read $ args !! 3
        e = read $ args !! 4
        resname = args !! 5
    zmatrix <- readZMatrix zmatr_name
    molecule <- readMolecule mol_name
    let moleculeM = setAtomWithOutOptimization n zmatrix molecule
        molecule' = fromMaybe (error "insertMolecule: setAtomWithOutOptimization returned nothing") moleculeM
    molecule'' <- setAtomWithOptimization s e zmatrix molecule'
    writeMolecule resname molecule''

geneticAlgorithm :: IO ()
geneticAlgorithm = do
    let evolution n sc (pop, all, best) = do
        -- | Производим один шаг эволюции.
        -- Если это первая популяция, то только
        -- считаем параметры @pop_x@
            pop' <- 
                if n == 1 
                then pop >>= flip computeLambda all
                else pop >>= selection >>= crossover >>= mutation >>= flip computeLambda all
            
            -- | Печатаем новую популяцию @pop_x'@ в @out_file@
            -- (определение @out_file@ смотри в модуле IO.hs)
            writeInFile ("Step " <> show n <> "\n") pop'
            
            -- | Дополняем множество всех особей @all_x@ 
            -- новой популяцией @pop_x'@
            let all' = all <> pop'
        
            -- | Ищем лучшую особь в новой популяции @pop_x'@ 
            -- и сравниваем её с лучшей среди всех в @all_x@.
            -- Тут же обновляем счетчик @sc@.
            let (best', sc') = if maximum pop' > best 
                then (maximum pop', 0)
                else (best, sc + 1) 
        
            -- | Проверяем условия остановки.
            -- В зависимости от результата возвращаем
            -- все особи @all_x'@ или продолжаем работу.
            if isStop sc' then return all'
            else evolution (n + 1) sc' (return pop', all', best')
    res <- evolution 1 0 (generatePopulation, [], tmpProtein)
    writeInFile "--------------\n--------------\n" []
    writeInFile "Result:\n" (reverse $ sortOn lambda res)
    return ()

main :: IO ()
main = do
    print "Hello World!"