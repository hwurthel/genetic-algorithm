-- Модуль содержит реализацию процессов эволюции
module Evolution 
    ( generatePopulation
    , crossover
    , mutation
    , selection
    , computeLambda
    , writeInProteinFile
    ) where
        
import Protein
import Config (eval_param, pop_size, 
               prob_cros, prob_mut,)
import Config (out_file, tmp_of, tmp_if,time_wait)

import Control.Monad (replicateM)
import System.Random (randomRIO)
import System.IO
import System.Directory (removeFile, doesFileExist)
import Control.Concurrent
import Data.List
import Data.Maybe (fromJust)

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
    v <- mapM selectAminoacid bros_variance
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
     (for_cros, other) <- makeParentsPair <$> selectProtein prob_cros ps
     crossed <- mapM crossover' for_cros
     return $ uncurry (<>) $ revMakeParentsPair (crossed, other)

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
    mutated <- mapM mutation' for_mutation
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
selection p = selection' 0 (sortPopulation p) 
    where 
        q = map (\n -> sum $ map eval [1..n]) [1..pop_size]
        selectNum r (x:xs) = 1 + (if r > x then selectNum r xs else 0)
        selection' n p
            | n == pop_size = return []
            | otherwise = 
                do
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
    ps' <- zip ps <$> (replicateM (length ps) $ randomRIO (0, 1 :: Double))
    let (p1, p2) = partition (\x -> snd x < prob) ps'
    return (map fst p1, map fst p2)

-- | НАЧАЛО. ЧТЕНИЕ И ВЫВОД.
computeLambda :: [Protein] -> [Protein] -> IO [Protein]
computeLambda p ps = do
    -- |Реализовать алгоритм выбора белков @p_a@ и @p_b@
    -- можно иначе -- сначала найти @p_b@, а потом искать @p_a@
    -- как дополнение @p_b@ до @p@. Такая реализация
    -- окажется быстрее. СДЕЛАТЬ ПОТОМ
    let p_a = [x | x <- p , x `notElem` ps]
        p_b = map (\x -> x {lambda = lambda (inPs x)}) (p \\ p_a)
              where inPs x = fromJust $ find (== x) ps 
    
    print "Current population\n"
    print $ zip (map variance p) (map lambda p)     
    print "\nWill be computed"
    print $ zip (map variance p_a) (map lambda p_a)
    print "\nWon't be computed"
    print $ zip (map variance p_b) (map lambda p_b)
    print "\nAll proteins"   
    print $ zip (map variance ps) (map lambda ps)
    print "----------------------------------------------\n"

    withFile tmp_of WriteMode (writeProtein p_a)
    wait True
    p_a' <- withFile tmp_if ReadMode (readProtein p_a)
    removeFile tmp_of
    removeFile tmp_if
    return (p_a' <> p_b) 
    where 
        writeProtein ps hdl = mapM_ (hPutStrLn hdl . protein) ps
        readProtein  ps hdl = mapM (\p' -> hGetLine hdl >>= return . (\x -> p' { lambda = Just x }). read) ps
        wait False = return () -- Рассмотреть функцию hWaitForInput
        wait True  = do 
            e <- doesFileExist tmp_if 
            if e then wait False
            else threadDelay time_wait >> wait True

writeInProteinFile :: String -> [Protein] -> IO [Protein]
writeInProteinFile msg ps = withFile out_file AppendMode write >> return ps
            where write = \hdl -> hPutStrLn hdl msg >> mapM_ (\p' -> hPrint hdl p') ps
-- | КОНЕЦ. ЧТЕНИЕ И ВЫВОД.