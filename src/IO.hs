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
import Data.List          (notElem, partition, intersectBy, find, (\\))
import Data.Maybe

import Protein

out_file  = "result"
tmp_of    = "tmp_out"
tmp_if    = "tmp_in"
time_wait = 5*10^6 -- микросекунды

computeLambda :: [Protein] -> [Protein] -> IO [Protein]
computeLambda ps p = do
    let p_a = [x | x <- p , (variance x) `notElem` (map variance ps)]
        p_b = p \\ p_a
        p_b' = if p_b == []
            then []
            else map (flip search ps) p_b
    print "P\n"
    print (map variance p)        
    print "P_A\n"
    print (map variance p_a)
    print "P_B\n"
    print (map variance p_b)
    print "P_B'\n"   
    print (map variance p_b')
    print "PS\n"   
    print (map variance ps)        
    print "Compute1"
    withFile tmp_of WriteMode (writeProtein p_a)
    print "Compute2"
    wait True
     
    p_a' <- withFile tmp_if ReadMode (readProtein p_a)
    removeFile tmp_of
    removeFile tmp_if
    return (p_a' <> p_b') 
    where 
        writeProtein ps = \hdl -> mapM_ (\p' -> hPutStrLn hdl $ protein p') ps
        readProtein  ps = \hdl -> mapM (\p' -> hGetLine hdl >>= return . (update p'). read) ps
        update p x = p { lambda = Just x }
        search x y = Protein {variance = variance z, lambda = lambda z, protein = protein z}
            where z = fromJust $ find (\y' -> variance y' == variance x) y

        wait False = return ()
        wait True  = do 
            e <- doesFileExist tmp_if 
            if e then wait False
            else threadDelay time_wait >> wait True

writeInFile :: String -> [Protein] -> IO [Protein]
writeInFile msg ps = withFile out_file AppendMode writeProtein >> return ps
            where writeProtein = \hdl -> hPutStrLn hdl msg >> mapM_ (\p' -> hPutStrLn hdl $ show p') ps