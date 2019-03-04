-- Модуль содержит описание структуры Protein
module Protein where

type Aminoacid = Char

-- Структура белка
data Protein = Protein { protein   :: [Aminoacid]
                       , variance  :: [Aminoacid]
                       , lambda    :: Double
                       } 

-- Полагаем, что белок является "лучше", если
-- значение lambda больше
instance Eq Protein where
    (==) a b = (==) (lambda a) (lambda b)

instance Ord Protein where
    compare a b = compare (lambda a) (lambda b)

instance Show Protein where
    show p = "\nLambda:\t" <> (show $ lambda p) <>
             "\nVariance:\t" <> (variance p) <>
             "\nProtein:\t" <> (protein p)

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
    lambda   = 750
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


-- Вставка изменяемых аминокислот в шаблонный белок
-- с целью получить вид полученного белка  
insertVariance :: [(Aminoacid, Int)] -> [Aminoacid]
insertVariance x = foldl insert (protein tmpProtein) x
    where insert p (a, n) = take (n-1) p <> [a] <> drop n p