-- Модуль содержит описание структуры Protein
module Protein where

import System.Random (randomRIO)

type Aminoacid = Char

-- Структура белка
data Protein = Protein { protein   :: [Aminoacid]
                       , variance  :: [Aminoacid]
                       , lambda    :: Maybe Double
                       } 

-- | Белки одинаковы, если последовательности
-- аминокислот одинакова @protein@. Но так как отличия 
-- только в @variance@, то сравниваем это поле.
instance Eq Protein where
    (==) a b = (==) (variance a) (variance b)

-- | Белок является "лучше", если
-- значение поля @lambda@ больше.
instance Ord Protein where
    compare a b = compare (lambda a) (lambda b)

instance Show Protein where
    show p = "Lambda:    \t" <> (show $ lambda p) <>
             "\nVariance:\t" <> (variance p)      <>
             "\nProtein: \t" <> (protein p)       <>
             "\n"

-- Шаблонный белок
tmpProtein :: Protein
tmpProtein = Protein { 
    variance = [],
    protein  = "MLMTVFSSAPELALLGSTFAQVDPSNLSVSDSLTYGQFNLVYNA" <>
               "FSFAIAAMFASALFFFSAQALVGQRYRLALLVSAIVVSIAGYHY" <>
               "FRIFNSWDAAYVLENGVYSLTSEKFNDAYRYVDWLLTVPLLLVE" <>
               "TVAVLTLPAKEARPLLIKLTVASVLMIATGYPGEISDDITTRII" <>
               "WGTVSTIPFAYILYVLWVELSRSLVRQPAAVQTLVRNMRWLLLL" <>
               "SWGVYPIAYLLPMLGVSGTSAAVGVQVGYTIADVLAKPVFGLLV" <>
               "FAIALVKTKADQESSEPHAAIGAAANKSGGSLIS",
    lambda   = Nothing
    }

    
-- Пары, определяющие, в каком месте белка и на какую аминокислоту
-- мы можем произвести замену
bros :: [([Aminoacid], Int)]
bros = [("DAPQGSKTLVNWM", 121),
        ("TCDGLNV"      , 125),
        ("LKW"          , 129),
        ("WM"           , 222),
        ("YA"           , 225),
        ("PIV"          , 226),
        ("DET"          , 253),
        ("AM"           , 256)]

bros_var = fst . unzip $ bros
bros_pos = snd . unzip $ bros

-- Вставка изменяемых аминокислот в шаблонный белок
-- с целью получить вид полученного белка  
insertVariance :: [(Aminoacid, Int)] -> [Aminoacid]
insertVariance x = foldl insert (protein tmpProtein) x
    where insert p (a, n) = take (n-1) p <> [a] <> drop n p

-- Взятие случайной аминокислоты из набора доступных
selectAminoacid :: [Aminoacid] -> IO Aminoacid
selectAminoacid x = do
    r <- randomRIO (0, length x - 1)
    return $ x !! r

-- Расчет вероятности того, что будет сгенерирована 
-- полная популяция при заданных числах генов и особей 
-- при указанных выше параметрах
getProbFulPop :: [Integer] -> [(Double, Integer)]
getProbFulPop ns = zip f ns
    where 
        f = map (product.(\n -> (map probFulPop l) <*> [n])) ns
        l = map (fromIntegral.length.fst) bros

-- Расчет вероятности того, что будет сгенерирована 
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