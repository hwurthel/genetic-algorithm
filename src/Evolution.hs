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
import Config ( evalParam
              , popSize
              , probCros
              , probMut
              , resultFile
              , computeLambdaOuf
              , computeLambdaInf
              , timeWait)

import System.Directory (removeFile, doesFileExist, renameFile)
import Control.Monad (replicateM)
import System.Random (randomRIO)
import Data.Maybe (fromJust)
import Control.Concurrent
import System.IO
import Data.List
import Data.Default

-- | НАЧАЛО. ГЕНЕРАЦИЯ ПОПУЛЯЦИИ
-- | Инициализатор популяции
generatePopulation :: IO [Protein]
generatePopulation =
    let generatePopulation' 0 = [] 
        generatePopulation' n = generateProtein : generatePopulation' (n-1)
    in  sequence $ generatePopulation' popSize

-- | Инициализатор особи
generateProtein :: IO Protein
generateProtein = do
    v <- mapM selectAminoacid brosVariance
    return $ def { variance = v, 
                   protein  = insertVariance $ zip v brosPosition, 
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
     (forCros, other) <- makeParentsPair <$> selectProtein probCros ps
     crossed <- mapM crossover' forCros
     return $ uncurry (<>) $ revMakeParentsPair (crossed, other)

crossover' :: (Protein, Protein) -> IO (Protein, Protein)
crossover' (p1, p2) = do
     let (v1, v2) = (variance p1, variance p2)
     (v1', v2') <- crossover'' (v1, v2) brosProbCros
     let p1' = def {variance = v1', protein = insertVariance $ zip v1' brosPosition, lambda = Nothing }
         p2' = def {variance = v2', protein = insertVariance $ zip v2' brosPosition, lambda = Nothing }
     return (p1', p2')
     where
        crossover'' ([], []) []  = return ([], [])
        crossover'' ((x:xs), (y:ys)) (a:as) = do
                r <- randomRIO (0, 1 :: Double)
                if r < a 
                 then pure ([y], [x]) <> crossover'' (xs, ys) as
                 else pure ([x], [y]) <> crossover'' (xs, ys) as
                    
-- | Образовываем родительские пары. Кому-то может не достаться особи. 
-- Такая особь отправляется во вторую группу.
makeParentsPair :: ([Protein], [Protein]) -> ([(Protein, Protein)], [Protein])
makeParentsPair ([], a) = ([], a)
makeParentsPair ([y], a) = ([], y:a)
makeParentsPair (x:y:s, a) = ([(x,y)], []) <> makeParentsPair (s, a)

-- | Операция, обратная makeParentsPair
revMakeParentsPair :: ([(Protein, Protein)], [Protein]) -> ([Protein], [Protein])
revMakeParentsPair ([]  , y) = ([], y)
revMakeParentsPair (x:xs, y) = ([fst x], []) <> ([snd x], []) <> revMakeParentsPair (xs, y)
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
    (forMutation, other) <- selectProtein probMut ps
    mutated <- mapM mutation' forMutation
    return $ mutated <> other

mutation' :: Protein -> IO Protein
mutation' p = do
    v' <- mutation'' (variance p) brosVariance brosProbMut
    return $ Protein {variance = v', protein = insertVariance $ zip v' brosPosition, lambda = Nothing }
    where
        mutation'' [] [] [] = return [] 
        mutation'' (x:xs) (y:ys) (a:as) = do
            r <- randomRIO (0, 1 :: Double)
            if r < a 
             then sequence [selectAminoacid $ delete x y] <> mutation'' xs ys as
             else pure [x] <> mutation'' xs ys as
-- КОНЕЦ. МУТАЦИИ 

-- | НАЧАЛО. СЕЛЕКЦИЯ
-- | Выбирается pop_size лучших особей.
selection :: [Protein] -> IO [Protein]
selection [] = return []
selection p = selection' 0 (sortPopulation p) 
    where 
        q = map (\n -> sum $ map eval [1..n]) [1..popSize]
        selectNum r (x:xs) = 1 + (if r > x then selectNum r xs else 0)
        selection' n p
            | n == popSize = return []
            | otherwise = do
                    r <- randomRIO (0, last q)
                    let x = p !! (selectNum r q - 1)
                    pure [x] <> selection' (n+1) p


-- | Функция ранжирования
eval :: Int -> Double
eval n = evalParam*(1-evalParam)^(n-1)

-- | Ранжирование популяции от лучшей (на первом месте) к худшей
sortPopulation :: [Protein] -> [Protein]
sortPopulation p = reverse $ sortOn lambda p 
-- | КОНЕЦ. СЕЛЕКЦИЯ

-- | Отбор хромосом для некоторого процесса. 
-- Функция возвращает пару, где на первом месте хромосомы, участвующие в некотором процессе,
-- на втором - все оставшиеся. Первый параметр -- вероятность попасть в первую группу.
selectProtein :: Double -> [Protein] ->  IO ([Protein], [Protein])
selectProtein prob ps = do
    ps' <- zip ps <$> replicateM (length ps) (randomRIO (0, 1 :: Double))
    let (p1, p2) = partition (\x -> snd x < prob) ps'
    return (map fst p1, map fst p2)

-- | НАЧАЛО. ЧТЕНИЕ И ВЫВОД.
computeLambda :: [Protein] -> [Protein] -> IO [Protein]
computeLambda pop all = do
    -- |Реализовать алгоритм выбора белков @p_a@ и @p_b@
    -- можно иначе -- сначала найти @p_b@, а потом искать @p_a@
    -- как дополнение @p_b@ до @p@. Такая реализация
    -- окажется быстрее. СДЕЛАТЬ ПОТОМ
    let pa = [x | x <- pop , x `notElem` all]
        pb = map (\x -> x {lambda = lambda (inPs x)}) (pop \\ pa)
         where inPs x = fromJust $ find (== x) all 
    
    print "Current population"
    print $ zip (map variance pop) (map lambda pop)     
    print "Will be computed"
    print $ zip (map variance pa) (map lambda pa)
    print "Won't be computed"
    print $ zip (map variance pb) (map lambda pb)
    print "All proteins"   
    print $ zip (map variance all) (map lambda all)
    print "----------------------------------------------"

    (tmpNameOuf, tmpHandleOuf) <- openTempFile "." "temp"
    mapM_ (hPutStrLn tmpHandleOuf . protein) pa
    hClose tmpHandleOuf
    renameFile tmpNameOuf computeLambdaOuf
    wait True
    handleInf <- openFile computeLambdaInf ReadMode
    pa' <- mapM (\p' -> hGetLine handleInf >>= return . (\x -> p' { lambda = Just x}) . read) pa
    hClose handleInf
    removeFile computeLambdaOuf
    removeFile computeLambdaInf

    return (pa' <> pb) 
    where 
        wait False = return ()
        wait True  = do 
            e <- doesFileExist computeLambdaInf 
            if e then wait False
            else threadDelay timeWait >> wait True

writeInProteinFile :: String -> [Protein] -> IO [Protein]
writeInProteinFile msg ps = withFile resultFile AppendMode write >> return ps
            where write hdl = hPutStrLn hdl msg >> mapM_ (hPrint hdl) ps
-- | КОНЕЦ. ЧТЕНИЕ И ВЫВОД.