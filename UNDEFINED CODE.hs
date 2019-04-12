-- | Расчет вероятности того, что будет сгенерирована 
-- полная популяция при заданных числах генов и особей 
-- при указанных выше параметрах
getProbFulPop :: [Integer] -> [Integer] -> [(Double, Integer)]
getProbFulPop ls ns = zip f ns
    where 
        f = map (product (\n -> map probFulPop ls <*> [n]) ns

-- | Расчет вероятности того, что будет сгенерирована 
-- полная популяция при заданных числах генов и особей
-- m - число генов, n - число особей
factorial :: Integer -> Integer
factorial 0 = 1
factorial n = n * factorial (n - 1)

binom :: Integer -> Integer -> Integer
binom m n = ceiling $ a / b
        where a = realToFrac $ factorial m
              b = realToFrac $ factorial n * factorial (m - n)

probFulPop :: Integer -> Integer -> Double
probFulPop m n
    | m > n = 0
    | otherwise = 1 - realToFrac (probFulPop' m n 1) / realToFrac m^n
    where
        probFulPop' 1 _ _ = 0
        probFulPop' m n s = ceiling $ compute
            where compute = realToFrac ((binom m (m - 1)) * ((m - 1)^n - probFulPop' (m - 1) n (s + 1))) / realToFrac (factorial s)


              
-- findRings :: Molecule -> [[Index]]
-- findRings molecule =
--     let index = Map.keys  $ getBonds molecule
--         bonds = Map.elems $ getBonds molecule
--         filterEnum xs x
--             | (head xs == x) && (length xs > 2) = [xs ++ [x]]
--             | length xs > 7 = []
--             | x `elem` xs = []
--             | otherwise = createEnum (xs ++ [x]) $ bonds !! (x - 1)
--         createEnum xs1 xs2 = concat $ (\x -> filterEnum xs1 x) <$> xs2
--     in concat $ (\x -> createEnum [x] $ bonds !! (x - 1)) <$> index