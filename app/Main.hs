module Main where

import StopCondition
import Evolution
import Protein
import InputOutput

import Data.List

main :: IO ()
main = do
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
            all' <- return $ all <> pop'

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