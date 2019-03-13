module Main where

import StopCondition
import Evolution
import Protein
import IO

import Data.List (sortOn)

main :: IO ()
main = do
    let evolution n sc (pop_x, all_x, best_x) = do
            -- | Производим один эволюционный шаг.
            -- Если это первая популяция, то только
            -- считаем параметры @pop_x@
            if n == 1 
                then pop_x' <- pop_x >>= flip computeLambda all_x
                else pop_x' <- pop_x >>= selection >>= crossover >>= mutation >>= flip computeLambda all_x

            -- | Печатаем новую популяцию @pop_x'@ в @out_file@
            -- (определение имени файла смотреть в модуле IO.hs)
            writeInFile ("Step " <> show n <> "\n") pop_x'
            
            -- | Дополняем множество всех особей @all_x@ 
            -- новой популяцией @pop_x'@
            all_x' <- return $ all_x <> pop_x'

            -- | Ищем лучшую особь в новой популяции @pop_x'@ 
            -- и сравниваем её с лучшей среди всех в @all_x@.
            -- Тут же обновляем счетчик @sc@.
            let (best_x', sc') = if maximum pop_x' > best_x 
                then (maximum pop_x', 0)
                else (best_x, sc + 1) 

            -- | Проверяем условия остановки.
            -- В зависимости от результата возвращаем
            -- все особи @all_x'@ или продолжаем работу.
            if isStop sc' then return all_x'
            else process (n + 1) sc' (return pop_x', all_x', best_x')
            
    res <- process 1 0 (generatePopulation, [], tmpProtein)
    writeInFile "--------------\n--------------\n" []
    writeInFile "Result:\n" (reverse $ sortOn lambda res)
    
    return ()