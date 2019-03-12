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
import Data.List          (notElem, partition, intersectBy)

import Protein

out_file  = "result"
tmp_of    = "tmp_out"
tmp_if    = "tmp_in"
time_wait = 5*10^6 -- микросекунды

computeLambda :: [Protein] -> [Protein] -> IO [Protein]
computeLambda p ps = do
    let (p_a, p_b)   = ( [x | x <- p , variance x `notElem` (map variance ps)]
                       , [x | x <- ps, variance x `elem`    (map variance p) ] )
    withFile tmp_of WriteMode (writeProtein p_a)
    wait True 
    p_a' <- withFile tmp_if ReadMode (readProtein p_a)
    removeFile tmp_of
    removeFile tmp_if
    return (p_a' <> p_b) 
    where 
        writeProtein ps = \hdl -> mapM_ (\p' -> hPutStrLn hdl $ protein p') ps
        readProtein  ps = \hdl -> mapM (\p' -> hGetLine hdl >>= return . (update p'). read) ps
        update p x = p { lambda = Just x }
        wait False = return ()
        wait True  = do 
            e <- doesFileExist tmp_if 
            if e then wait False
            else threadDelay time_wait >> wait True

writeInFile :: String -> [Protein] -> IO [Protein]
writeInFile msg ps = withFile out_file AppendMode writeProtein >> return ps
            where writeProtein = \hdl -> hPutStrLn hdl msg >> mapM_ (\p' -> hPutStrLn hdl $ show p') ps
