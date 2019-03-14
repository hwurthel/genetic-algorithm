-- Модуль содержит реализацию процессов эволюции
module Evolution 
    (
      generatePopulation
    , crossover
    , mutation
    , selection
    ) where
        
import Protein
import Config (eval_param, pop_size, prob_cros, 
        prob_cros_gene, prob_mut, prob_mut_gene)

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
    v <- sequence $ fmap selectAminoacid bros_var
    return $ Protein { variance = v, 
                       protein  = insertVariance $ zip v bros_pos, 
                       lambda   = Nothing}
-- | КОНЕЦ. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ

-- | НАЧАЛО. КРОССИНГОВЕР
-- | Кроссинговер между особями на множестве особей
crossover :: [Protein] -> IO [Protein]
crossover [] = return []
crossover ps = do
     (ps'_f, ps'_s) <- selectProtein prob_cros ps >>= return . makeParentsPair
     ps'_cros       <- sequence $ map crossover' ps'_f
     return $ (\(x,y) -> x <> y) (revMakeParentsPair (ps'_cros, ps'_s))

-- | Кроссинговер между парой особей
crossover' :: (Protein, Protein) -> IO (Protein, Protein)
crossover' (p1, p2) = do
     let (v1, v2) = (variance p1, variance p2)
     (v1', v2') <- crossover'' (v1, v2) 0
     let p1' = Protein {variance = v1', protein = insertVariance $ zip v1' bros_pos, lambda = Nothing }
         p2' = Protein {variance = v2', protein = insertVariance $ zip v2' bros_pos, lambda = Nothing }
     return (p1', p2')

-- | Кроссинговер между парой особей. Вспомогательная функция
crossover'' :: ([Aminoacid], [Aminoacid]) -> Int -> IO ([Aminoacid], [Aminoacid]) 
crossover'' (v1, v2) n
    | n == length bros = return (v1, v2)
    | otherwise = do
        a <- randomRIO (0, 1 :: Double)
        if a > prob_cros_gene then crossover'' (v1, v2) (n + 1)
        else crossover'' (v1', v2') (n + 1)
            where v1' = take n v1 <> [v2 !! n] <> drop (n + 1) v1
                  v2' = take n v2 <> [v1 !! n] <> drop (n + 1) v2
                    
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
-- Каждый ген в особь подвергается изменению аминокислоты
-- из набора изменяемых, исключая текущую аминоксилоты,  
-- с вероятностью @prob_mut_gene@.
mutation :: [Protein] -> IO [Protein]
mutation [] = return []
mutation ps = do
    (ps'_f, ps'_s) <- selectProtein prob_mut ps
    ps'_mut        <- sequence $ map mutation' ps'_f
    return (ps'_mut <> ps'_s) 

mutation' :: Protein -> IO Protein
mutation' p = do
    v' <- sequence $ map (uncurry mutation'') (zip (variance p) bros_var)
    let p' = Protein {variance = v', protein = insertVariance $ zip v' bros_pos, lambda = Nothing }
    return p'

mutation'' :: Aminoacid -> [Aminoacid] -> IO Aminoacid
mutation'' a as = do
    r <- randomRIO (0, 1 :: Double)
    if r > prob_mut_gene then return a
    else selectAminoacid (delete a as) 
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
                return [x] <> selection' (n+1) p


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

