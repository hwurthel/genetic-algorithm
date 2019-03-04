-- Модуль взаимодействия программы
-- с внешними программами
module IO where

import System.IO
import Protein

outFileName = "output"
inFileName  = "input"

computeLambda :: [Protein] -> IO [Protein]
computeLambda p = do
    hOut <- openFile outFileName WriteMode
    mapM_ (\p' -> hPutStrLn hOut $ protein p') p
    hClose hOut
    -- Ожидаем вычисленния результата другим модулем
    -- ...
    hIn <- openFile inFileName ReadMode
    res <- mapM (\p' -> hGetLine hIn >>= readIO >>= return . update p') p
    hClose hIn
    return res
    where update p x = p { lambda = x }