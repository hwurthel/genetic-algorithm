module Evolution where

import System.Random
import Data.List
import Protein

a         = 0.05  -- Параметр для функции оценки
pop_size  = 10    -- Размер популяции
prob_cros = 0.3   -- Вероятность кроссинговера
prob_mut  = 0.2   -- Вероятность мутации

tl = 551.0 -- Искомый парамтер lambda, нм

-- НАЧАЛО. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ
-- Инициализатор популяции
generatePopulation :: IO [Protein]
generatePopulation =
    let generatePopulation' 0 = [] 
        generatePopulation' n = generateProtein : generatePopulation' (n-1)
    in  sequence $ generatePopulation' pop_size

-- Инициализатор особи
generateProtein :: IO Protein
generateProtein = do
    m' <- sequence $ fmap selectAminoacid bros
    return $ Protein { variance = fmap fst m', 
                       protein  = insertVariance m', 
                       lambda   = -1}
-- КОНЕЦ. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ

-- НАЧАЛО. КРОССИНГОВЕР
-- Кроссинговер между особями на множестве особей
crossingover :: [Protein] -> IO [Protein]
crossingover p = do
     p'       <- selectParents p
     pop_pair <- return   $ makeParentsPair p'
     pcf'     <- sequence $ map crossingover' (fst pop_pair)
     pcs'     <- return   $ snd pop_pair
     return $ (\(x,y) -> x <> y) (revMakeParentsPair (pcf',pcs') )

-- Кроссинговер между парой особей
-- Здесь надо подумать над тем, как можно реализовывать кроссинговер.
-- На данный момент алгоритм следующий: происходит обмен
-- аминокислотами с вероятностью 0.3.
-- Возможно, стоит реализовать следующую возможность:
-- после генерации потомка проверять его параметр lambda.
-- Если параметр стал хуже, то повторять кроссинговер.
-- Если параметр стал не хуже, то отправлять в следующую популяцию. 
crossingover' :: (Protein, Protein) -> IO (Protein, Protein)
crossingover' (x, y) = do
     (v1 , v2 ) <- return (variance x, variance y)
     (v1', v2') <- crossingover'' (v1, v2) 0
     let x' = x {variance = v1', lambda = -1 } 
         y' = y {variance = v2', lambda = -1 }
     return (x', y')
     where crossingover'' (v1, v2) n
            | n == length bros = return (v1, v2)
            | otherwise = do
                a <- randomRIO (0, 1 :: Double)
                if a > 0.3 then crossingover'' (v1, v2) (n + 1)
                else do
                    let v1' = take n v1 <> [v2 !! n] <> drop (n + 1) v1
                        v2' = take n v2 <> [v1 !! n] <> drop (n + 1) v2
                    crossingover'' (v1', v2') (n + 1)
                    

-- Отбор случайным образом родительских хромосом, участвующих в кроссинговере 
-- На первом месте в паре -- будущие родители, на втором -- все оставшиеся
selectParents :: [Protein] -> IO ([Protein], [Protein])
selectParents p = selectParents' 0 p where
    selectParents' n p
         | n == pop_size = return ([],[])
         | otherwise = do
            a <- randomRIO (0, 1 :: Double)
            if a < prob_cros 
                then return ([p !! n], []) <> selectParents' (n+1) p 
                else return ([], [p !! n]) <> selectParents' (n+1) p

-- Образовываем родительские пары. Кому-то может не достаться особи. Такая особь отправляется во вторую группу.
makeParentsPair :: ([Protein], [Protein]) -> ([(Protein, Protein)], [Protein])
makeParentsPair ([], a) = ([], a)
makeParentsPair (y:[], a) = ([], y:a)
makeParentsPair (x:y:s, a) = ([(x,y)], []) <> makeParentsPair (s, a)

-- Операция, обратная makeParentsPair
revMakeParentsPair :: ([(Protein, Protein)], [Protein]) -> ([Protein], [Protein])
revMakeParentsPair ([]  , y) = ([], y)
revMakeParentsPair ((x:xs), y) = ([fst x], []) <> ([snd x], []) <> revMakeParentsPair (xs, y)
-- КОНЕЦ. КРОССИНГОВЕР

-- НАЧАЛО. МУТАЦИИ 
-- КОНЕЦ. МУТАЦИИ 

-- НАЧАЛО. СЕЛЕКЦИЯ
-- Целевая функция.
-- Здесь надо понимать, что явлется критерием того, что белок является "хорошм"
-- Допустим, что нам по условию задачи надо найти такой белок, который
-- поглощаем э/м волны на частоте tl нм
targetFun :: Protein -> Double
targetFun p = abs $ tl - lambda p

-- Ранжирование популяции от лучшей (на первом месте) к худшей
sortPopulation :: [Protein] -> [Protein]
sortPopulation p = sortOn targetFun p
-- КОНЕЦ. СЕЛЕКЦИЯ

-- Взятие случайной аминокислоты из набора доступных
selectAminoacid :: ([Aminoacid], Int) -> IO (Aminoacid, Int)
selectAminoacid (a, n)= do
    r <- randomRIO (0, length a - 1)
    return (a !! r, n)