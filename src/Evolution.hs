-- Модуль содержит реализацию процессов эволюции
module Evolution 
    (
      generatePopulation
    , crossover
    , mutation
    , selection
    ) where
        
import Protein
import Config (eval_param, pop_size, 
               prob_cros, prob_mut,)

import System.Random
import Data.List

-- | НАЧАЛО. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ
-- | Инициализатор популяции
generatePopulation :: IO [Protein]
generatePopulation =
    let generatePopulation' 0 = [] 
        generatePopulation' n = generateProtein : generatePopulation' (n-1)
    in  sequence $ generatePopulation' pop_size

-- | Инициализатор особи
generateProtein :: IO Protein
generateProtein = do
    v <- sequence $ fmap selectAminoacid bros_variance
    return $ Protein { variance = v, 
                       protein  = insertVariance $ zip v bros_position, 
                       lambda   = Nothing}
-- | КОНЕЦ. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ

-- | НАЧАЛО. КРОССИНГОВЕР
-- | На данный момент алгоритм следующий: вероятность того,
-- что особь подвергнется кроссоверу, определяется @prob_cros@.
-- Каждая аминокислота подвергается кроссоверу  
-- с соответствующей вероятностью, определяемой в @bros_prob_cros@.
crossover :: [Protein] -> IO [Protein]
crossover [] = return []
crossover ps = do
     (for_cros, other) <- selectProtein prob_cros ps >>= return . makeParentsPair
     crossed <- sequence $ map crossover' for_cros
     return $ (\(x,y) -> x <> y) $ revMakeParentsPair (crossed, other)

crossover' :: (Protein, Protein) -> IO (Protein, Protein)
crossover' (p1, p2) = do
     let (v1, v2) = (variance p1, variance p2)
     (v1', v2') <- crossover'' (v1, v2) bros_prob_cros
     let p1' = Protein {variance = v1', protein = insertVariance $ zip v1' bros_position, lambda = Nothing }
         p2' = Protein {variance = v2', protein = insertVariance $ zip v2' bros_position, lambda = Nothing }
     return (p1', p2')
     where
        crossover'' ([], []) []  = return ([], [])
        crossover'' ((x:xs), (y:ys)) (a:as) = do
                r <- randomRIO (0, 1 :: Double)
                if r < a then pure ([y], [x]) <> crossover'' (xs, ys) as
                         else pure ([x], [y]) <> crossover'' (xs, ys) as
                    
-- | Образовываем родительские пары. Кому-то может не достаться особи. 
-- Такая особь отправляется во вторую группу.
makeParentsPair :: ([Protein], [Protein]) -> ([(Protein, Protein)], [Protein])
makeParentsPair ([], a) = ([], a)
makeParentsPair (y:[], a) = ([], y:a)
makeParentsPair (x:y:s, a) = ([(x,y)], []) <> makeParentsPair (s, a)

-- | Операция, обратная makeParentsPair
revMakeParentsPair :: ([(Protein, Protein)], [Protein]) -> ([Protein], [Protein])
revMakeParentsPair ([]  , y) = ([], y)
revMakeParentsPair ((x:xs), y) = ([fst x], []) <> ([snd x], []) <> revMakeParentsPair (xs, y)
-- | КОНЕЦ. КРОССИНГОВЕР

-- | НАЧАЛО. МУТАЦИИ
-- | На данный момент алгоритм следующий: вероятность того,
-- что особь подвергнется мутации, определяется @prob_mut@.
-- Каждая аминокислота подвергается изменению
-- из набора изменяемых, исключая текущую аминокислоту,  
-- с соответствующей вероятностью, определяемой в @bros_prob_mut@.
mutation :: [Protein] -> IO [Protein]
mutation [] = return []
mutation ps = do
    (for_mutation, other) <- selectProtein prob_mut ps
    mutated <- sequence $ map mutation' for_mutation
    return $ mutated <> other

mutation' :: Protein -> IO Protein
mutation' p = do
    v' <- mutation'' (variance p) bros_variance bros_prob_mut
    return $ Protein {variance = v', protein = insertVariance $ zip v' bros_position, lambda = Nothing }
    where
        mutation'' [] [] [] = return [] 
        mutation'' (x:xs) (y:ys) (a:as) = do
            r <- randomRIO (0, 1 :: Double)
            if r < a then sequence [selectAminoacid $ delete x y] <> mutation'' xs ys as
                     else pure [x] <> mutation'' xs ys as
-- КОНЕЦ. МУТАЦИИ 

-- | НАЧАЛО. СЕЛЕКЦИЯ
-- | Выбирается pop_size лучших особей.
selection :: [Protein] -> IO [Protein]
selection [] = return []
selection p = do
    selection' 0 (sortPopulation p) 
    where 
        q = map (\n -> sum $ map eval [1..n]) [1..pop_size]
        selectNum r (x:xs) = 1 + (if r > x then selectNum r xs else 0)
        selection' n p
            | n == pop_size = return []
            | otherwise = do
                r <- randomRIO (0, last q)
                let i = selectNum r q
                    x = p !! (i - 1)
                pure [x] <> selection' (n+1) p


-- | Функция ранжирования
eval :: Int -> Double
eval n = eval_param*(1-eval_param)^(n-1)

-- | Ранжирование популяции от лучшей (на первом месте) к худшей
sortPopulation :: [Protein] -> [Protein]
sortPopulation p = reverse $ sortOn lambda p 
-- | КОНЕЦ. СЕЛЕКЦИЯ

-- | Отбор хромосом для некоторого процесса. 
-- Функция возвращает пару, где на первом месте хромосомы, участвующие в некотором процессе,
-- на втором - все оставшиеся. Первый параметр -- вероятность попасть в первую группу.
selectProtein :: Double -> [Protein] ->  IO ([Protein], [Protein])
selectProtein prob ps = do
    ps' <- (sequence $ take (length ps) $ repeat $ randomRIO (0, 1 :: Double)) >>= return . zip ps
    let (p1, p2) = partition (\x -> snd x < prob) ps'
    return (map fst p1, map fst p2)

