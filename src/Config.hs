module Config where

import System.IO.Unsafe
import Data.List.Split
import Data.List

type Config = [(String, String)]

config_file = "config"

readConfig :: Config
readConfig = [(key, val) | 
                line <- filter (\x -> (not $ "--" `isPrefixOf` x) && x `notElem` [""]) $ lines s,
                let (key : _ : val : _) = words line]
    where s = unsafePerformIO $ readFile config_file 

get :: String -> (String -> a) -> a
get key f = case lookup key readConfig of
    Nothing -> error ("get: not found: " ++ key)
    Just x  -> f x

-- | Параметр @protein@ шаблонного белка
tmp_protein :: String
tmp_protein = get "tmp_protein" id

-- | Параметр @variance@ шаблонного белка
tmp_variance :: String
tmp_variance =  get "tmp_variance" id

-- | Параметр @lambda@ шаблонного белка
tmp_lambda :: Maybe Double
tmp_lambda = get "tmp_lambda" read

-- | Набор заменяемых аминокислот
bros_list :: [([Char], Int, Double, Double)]
bros_list = get "bros" readBros
    where readBros s = [(id a, read b, read c, read d) 
                        | p <- splitOn "_" s, 
                        let (a:b:c:d:[]) = splitOn ":" p]

-- | Выходной файл
out_file :: String
out_file = get "out_file" id 

-- | Входной файл для обсчета параметров белка
tmp_of :: String
tmp_of = get "tmp_of" id 

-- | Входной файл, где находятся обсчитанные параметры
tmp_if :: String
tmp_if = get "tmp_if" id

-- | Время задержки перед проверка файла @tmp_if@, в мкс
time_wait :: Int
time_wait = get "time_wait" read

-- | Параметр для функции оценки
eval_param :: Double
eval_param = get "eval_param" read

-- | Размер популяции
pop_size :: Int
pop_size = get "pop_size" read

-- | Вероятность того, что хромосома будет участвовать в кроссинговере
prob_cros :: Double 
prob_cros = get "prob_cros" read

-- | Вероятность того, что хромосома будет участвовать в мутации
prob_mut :: Double
prob_mut = get "prob_mut" read