module Evolution where

import System.Random
import Data.List
import Protein

a  = 0.05  -- Параметр для функции оценки
ps = 10    -- Размер популяции
pc = 0.3   -- Вероятность кроссинговера
pm = 0.2   -- Вероятность мутации

-- НАЧАЛО. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ
-- Инициализатор популяции
generatePopulation :: IO [Protein]
generatePopulation =
    let generatePopulation' 0 = [] 
        generatePopulation' n = generateProtein : generatePopulation' (n-1)
    in  sequence $ generatePopulation' ps

-- Инициализатор особи
generateProtein :: IO Protein
generateProtein = do
    m' <- sequence $ fmap selectAminoacid bros
    return $ Protein { variance = fmap fst m', 
                       protein  = insertVariance m', 
                       lambda   = -1}
-- КОНЕЦ. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ

-- НАЧАЛО. КРОССИНГОВЕР 
-- КОНЕЦ. КРОССИНГОВЕР

-- НАЧАЛО. МУТАЦИИ 
-- КОНЕЦ. МУТАЦИИ 

-- НАЧАЛО. СЕЛЕКЦИЯ
-- Целевая функция.
-- Здесь надо понимать, что явлется критерием того, что белок является "хорошм"
-- Допустим, что нам по условию задачи надо найти такой белок, который
-- поглощаем э/м волны на частоте 760 нм
targetFun :: Protein -> Double
targetFun p = abs $ 760 - lambda p

-- Ранжирование популяции от лучшей (на первом месте) к худшей
sortPopulation :: [Protein] -> [Protein]
sortPopulation p = sortOn targetFun p
-- КОНЕЦ. СЕЛЕКЦИЯ

-- Взятие случайной аминокислоты из набора доступных
selectAminoacid :: ([Aminoacid], Int) -> IO (Aminoacid, Int)
selectAminoacid (a, n)= do
    r <- randomRIO (0, length a - 1)
    return (a !! r, n)