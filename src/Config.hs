module Config where

import System.IO.Unsafe
import Data.List.Split
import Data.List

type Config = [(String, String)]

configFile = "config"

readConfig :: Config
readConfig = 
    [(key, val) | 
    line <- filter (\x -> (not $ "--" `isPrefixOf` x) && x /= "" ) $ lines s,
    let (key, val) = (\x -> (head x, unwords $ drop 2 x)) $ words line]
    where 
        s = unsafePerformIO $ readFile configFile 
    -- readConfig = 
    --     [(key, val) | 
    --     line <- filter (\x -> (not $ "--" `isPrefixOf` x) && x /= "" ) $ lines s,
    --     let (key, val) = (\x -> (head x, unwords $ drop 2 x)) $ words line]
    --     where s = unsafePerformIO $ readFile configFile

get :: String -> (String -> a) -> a
get key f = case lookup key readConfig of
    Nothing -> error ("get: not found: " ++ key)
    Just x  -> f x

-- | Параметр @protein@ шаблонного белка
tmpProtein :: String
tmpProtein = get "tmp_protein" id

-- | Параметр @variance@ шаблонного белка
tmpVariance :: String
tmpVariance =  get "tmp_variance" id

-- | Параметр @lambda@ шаблонного белка
tmpLambda :: Maybe Double
tmpLambda = get "tmp_lambda" read

-- | Набор заменяемых аминокислот
brosList :: [(String, Int, Double, Double)]
brosList = get "bros" readBros
    where readBros s = [(id a, read b, read c, read d) 
                        | p <- splitOn "," s, 
                        let [a,b,c,d] = words p]

-- | Выходной файл, содержащий все белки
resultFile :: String
resultFile = get "result_file" id 

-- | Выходной файл (для обсчета параметров белка)
computeLambdaOuf :: String
computeLambdaOuf = get "compute_lambda_ouf" id 

-- | Входной файл, где находятся обсчитанные параметры
computeLambdaInf :: String
computeLambdaInf = get "compute_lambda_inf" id

-- | Время задержки перед проверка файла @tmp_if@, в мкс
timeWait :: Int
timeWait = get "time_wait" read

-- | Параметр для функции оценки
evalParam :: Double
evalParam = get "eval_param" read

-- | Размер популяции
popSize :: Int
popSize = get "pop_size" read

-- | Вероятность того, что хромосома будет участвовать в кроссинговере
probCros :: Double 
probCros = get "prob_cros" read

-- | Вероятность того, что хромосома будет участвовать в мутации
probMut :: Double
probMut = get "prob_mut" read