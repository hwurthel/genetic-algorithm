module Config where

import System.IO.Unsafe

type Config = [(String, String)]

config_file = "config"

readConfig :: Config
readConfig = [(key, val) | line <- filter (`notElem` [""]) $ lines s, let (key : _ : val : _) = words line]
    where s = unsafePerformIO $ readFile config_file 

get :: String -> (String -> a) -> a
get key f = case lookup key readConfig of
    Nothing -> error ("get: not found: " ++ key)
    Just x  -> f x

out_file, tmp_of, tmp_if :: String
out_file  = get "out_file" id                   -- Имя выходного файла
tmp_of    = get "tmp_of" id                     -- Имя выходной файла для обсчета параметров белка
tmp_if    = get "tmp_if" id                     -- Имя входного файла, где находятся обсчитанные параметры

time_wait :: Int
time_wait = get "time_wait" read                -- Время задержки перед проверка файла @tmp_if@

eval_param                :: Double  
pop_size                  :: Int
prob_cros, prob_cros_gene :: Double 
prob_mut, prob_mut_gene   :: Double 
eval_param      = get "eval_param" read         -- Параметр для функции оценки
pop_size        = get "pop_size" read           -- Размер популяции
prob_cros       = get "prob_cros" read          -- Вероятность того, что хромосома будет участвовать в кроссинговере
prob_cros_gene  = get "prob_cros_gene" read     -- Вероятность того, что ген в хромосоме подвергнется кроссинговеру
prob_mut        = get "prob_mut" read           -- Вероятность того, что хромосома будет участвовать в мутации
prob_mut_gene   = get "prob_mut_gene" read      -- Вероятность того, что ген в хромосоме подвергнется мутации