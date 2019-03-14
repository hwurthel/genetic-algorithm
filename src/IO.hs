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
import Data.List          (notElem, partition, intersectBy, 
                           find, (\\))
import Data.Maybe         (fromJust)

import Protein

out_file  = "result"
tmp_of    = "tmp_out"
tmp_if    = "tmp_in"
time_wait = 5*10^6 -- микросекунды

computeLambda :: [Protein] -> [Protein] -> IO [Protein]
computeLambda p ps = do
    -- |Реализовать алгоритм выбора белков @p_a@ и @p_b@
    -- можно иначе -- сначала найти @p_b@, а потом искать @p_a@
    -- как дополнение @p_b@ до @p@. Такая реализация
    -- окажется быстрее. СДЕЛАТЬ ПОТОМ
    let p_a = [x | x <- p , x `notElem` ps]
        p_b = map (\x -> x {lambda = lambda (inPs x)}) (p \\ p_a)
              where inPs = \x -> fromJust $ find (== x) ps 
    
    putStrLn "Current population"
    putStrLn $ show  (zip (map variance p) (map lambda p))        
    putStrLn "Will be computed"
    putStrLn $ show (zip (map variance p_a) (map lambda p_a))
    putStrLn "Won't be computed"
    putStrLn $ show (zip (map variance p_b) (map lambda p_b))
    putStrLn "All proteins"   
    putStrLn $ show (zip (map variance ps) (map lambda ps))
    putStrLn ""

    withFile tmp_of WriteMode (writeProtein p_a)
    wait True
    p_a' <- withFile tmp_if ReadMode (readProtein p_a)
    removeFile tmp_of
    removeFile tmp_if
    return (p_a' <> p_b) 
    where 
        writeProtein ps = \hdl -> mapM_ (\p' -> hPutStrLn hdl $ protein p') ps
        readProtein  ps = \hdl -> mapM (\p' -> hGetLine hdl >>= 
                                  return . (\x -> p' { lambda = Just x }). read) ps
        wait False = return ()
        wait True  = do 
            e <- doesFileExist tmp_if 
            if e then wait False
            else threadDelay time_wait >> wait True

writeInFile :: String -> [Protein] -> IO [Protein]
writeInFile msg ps = withFile out_file AppendMode writeProtein >> return ps
            where writeProtein = \hdl -> hPutStrLn hdl msg >> 
                                      mapM_ (\p' -> hPutStrLn hdl $ show p') ps