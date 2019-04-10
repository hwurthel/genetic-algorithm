module Main where

import StopCondition
import Evolution
import Protein
import InsertMolecule

import System.Environment
import Data.List (sortOn)
import Data.Maybe (fromMaybe)

insertMolecule :: IO ()
insertMolecule = do
    args <- getArgs
    let zmatr_name = args !! 1
        mol_name   = args !! 2
        n = read $ args !! 3
        s = read $ args !! 4
        e = read $ args !! 5
        resname = args !! 6
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
            writeInProteinFile ("Step " <> show n <> "\n") pop'
            
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
    writeInProteinFile "--------------\n--------------\n" []
    writeInProteinFile "Result:\n" (reverse $ sortOn lambda res)
    return ()

main :: IO ()
main = do
    args <- getArgs
    let programm = args !! 0
    case programm of
        "insertMolecule"   -> insertMolecule
        "geneticAlgorithm" -> geneticAlgorithm
    return ()