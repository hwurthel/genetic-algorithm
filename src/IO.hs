-- Модуль взаимодействия программы
-- с внешними программами
module IO where
    
import System.IO
import System.Directory   (removeFile, doesFileExist)
import Control.Concurrent (threadDelay)

import Protein

outFileName = "output"
inFileName  = "input"
timeWait    = 5*10^6 -- микросекунды

computeLambda :: [Protein] -> IO [Protein]
computeLambda p = do
    hOut <- openFile outFileName WriteMode
    mapM_ (\p' -> hPutStrLn hOut $ protein p') p
    hClose hOut
    hIn <- wait True (openFile inFileName ReadMode)
    res <- mapM (\p' -> hGetLine hIn >>= readIO >>= return . update p') p
    hClose hIn
    removeFile outFileName
    removeFile inFileName
    return res
    where 
        update p x = p { lambda = Just x }
        wait False f = f
        wait True  f = do
            e <- doesFileExist inFileName 
            if e then wait False f
            else do threadDelay timeWait
                    wait True f