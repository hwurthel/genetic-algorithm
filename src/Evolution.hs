module Evolution where

import System.Random (randomRIO)
import Data.List     (sortOn)
import Protein

a         = 0.05  -- Параметр для функции оценки
pop_size  = 100    -- Размер популяции
prob_cros = 0.3   -- Вероятность кроссинговера
prob_mut  = 0.2   -- Вероятность мутации

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
    m' <- sequence $ fmap selectAminoacid bros_var
    return $ Protein { variance = m', 
                       protein  = insertVariance $ zip m' bros_pos, 
                       lambda   = Nothing}
-- КОНЕЦ. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ

-- НАЧАЛО. КРОССИНГОВЕР
-- Кроссинговер между особями на множестве особей
crossingover :: [Protein] -> IO [Protein]
crossingover p = do
     p'       <- selectParents p
     p_pair   <- return   $ makeParentsPair p'
     pcf'     <- sequence $ map crossingover' (fst p_pair)
     pcs'     <- return   $ snd p_pair
     return $ (\(x,y) -> x <> y) (revMakeParentsPair (pcf',pcs') )

-- Кроссинговер между парой особей
crossingover' :: (Protein, Protein) -> IO (Protein, Protein)
crossingover' (x, y) = do
     (v1 , v2 ) <- return (variance x, variance y)
     (v1', v2') <- crossingover'' (v1, v2) 0
     let x' = x {variance = v1', protein = insertVariance $ zip v1' bros_pos, lambda = Nothing }
         y' = y {variance = v2', protein = insertVariance $ zip v2' bros_pos, lambda = Nothing }
     return (x', y')

-- Кроссинговер между парой особей. Вспомогательная функция
crossingover'' :: ([Aminoacid], [Aminoacid]) -> Int -> IO ([Aminoacid], [Aminoacid]) 
crossingover'' (v1, v2) n
    | n == length bros = return (v1, v2)
    | otherwise = do
        a <- randomRIO (0, 1 :: Double)
        if a > 0.3 then crossingover'' (v1, v2) (n + 1)
        else do
            let v1' = take n v1 <> [v2 !! n] <> drop (n + 1) v1
                v2' = take n v2 <> [v1 !! n] <> drop (n + 1) v2
            crossingover'' (v1', v2') (n + 1)
                    

-- Отбор случайным образом родительских хромосом, участвующих в кроссинговере 
-- На первом месте в паре - будущие родители, на втором - все оставшиеся
selectParents :: [Protein] -> IO ([Protein], [Protein])
selectParents p = selectParents' 0 p where
    selectParents' n p
         | n == pop_size = return ([],[])
         | otherwise = do
            a <- randomRIO (0, 1 :: Double)
            if a < prob_cros 
                then return ([p !! n], []) <> selectParents' (n+1) p 
                else return ([], [p !! n]) <> selectParents' (n+1) p

-- Образовываем родительские пары. Кому-то может не достаться особи. 
-- Такая особь отправляется во вторую группу.
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
-- На данный момент алгоритм следующий: мутация в одной аминокислоты
-- (из набора изменяемых, исключая текущую аминоксилоты) происходит с 
-- вероятностью 0.3.
-- Возможно, стоит реализовать следующую возможность:
-- после мутации проверять параметр lambda.
-- Если параметр стал хуже, то повторять мутацию.
-- Если параметр стал не хуже, то отправлять в следующую популяцию. 
mutation :: [Protein] -> IO [Protein]
mutation p = return p

mutation' :: Protein -> IO Protein
mutation' p = return p

mutation'' :: [Aminoacid] -> Int -> IO [Aminoacid]
mutation'' p _ = return p

-- КОНЕЦ. МУТАЦИИ 

-- НАЧАЛО. СЕЛЕКЦИЯ
-- Выбирается pop_size лучших особей.
selection :: [Protein] -> IO [Protein]
selection p = selection' 0 (sortPopulation p) 
    where selection' n p
            | n == pop_size = return []
            | otherwise = do
                let q = map cumul [1..pop_size]
                a <- randomRIO (0, last q)
                let i = select_num a q
                    x = p !! (i - 1)
                return [x] <> selection' (n+1) p
                where
                    cumul n = sum $ map eval [1..n]
                    select_num a (x:xs) = 1 + (if a > x then select_num a xs else 0)

-- Функция ранжирования
eval :: Int -> Double
eval n = a*(1-a)^(n-1)

-- Ранжирование популяции от лучшей (на первом месте) к худшей
sortPopulation :: [Protein] -> [Protein]
sortPopulation p = reverse $ sortOn lambda p 
-- КОНЕЦ. СЕЛЕКЦИЯ