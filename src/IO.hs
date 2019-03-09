-- Модуль взаимодействия программы
-- с внешними программами
module IO 
    (
      computeLambda
    , writeInFile
    ) where
    
import System.IO
import System.Directory   (removeFile, doesFileExist)
import Control.Concurrent (threadDelay)

import Protein

out_file = "oufile"
tmp_of = "tmp_out"
tmp_if  = "tmp_in"
time_wait    = 5*10^6 -- микросекунды

computeLambda :: [Protein] -> IO [Protein]
computeLambda ps = do
    withFile tmp_of WriteMode writeProtein
    wait True 
    res <- withFile tmp_if ReadMode readProtein
    removeFile tmp_of
    removeFile tmp_if
    return res
    where 
        writeProtein = \hdl -> mapM_ (\p' -> hPutStrLn hdl $ protein p') ps
        readProtein = \hdl -> mapM (\p' -> hGetLine hdl >>= return . (update p'). read) ps
        update p x = p { lambda = Just x }
        wait False = return ()
        wait True  = do 
            e <- doesFileExist tmp_if 
            if e then wait False
            else threadDelay time_wait >> wait True

writeInFile :: String -> [Protein] -> IO [Protein]
writeInFile msg ps = withFile out_file AppendMode writeProtein >> return ps
            where writeProtein = \hdl -> hPutStrLn hdl msg >> mapM_ (\p' -> hPutStrLn hdl $ show p') ps