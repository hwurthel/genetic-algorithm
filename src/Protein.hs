-- | Модуль содержит описание структуры Protein
module Protein where

import Data.List
import Utils
import Config (tmp_protein, tmp_lambda, tmp_variance, bros_list)
import System.Random

type Aminoacid = Char
data Protein = Protein { protein   :: [Aminoacid]
                       , variance  :: [Aminoacid]
                       , lambda    :: Maybe Double
                       } 

-- | Белки одинаковы, если последовательности
-- аминокислот @protein@ одинакова. Но так как отличия 
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

-- | Шаблонный белок
tmpProtein :: Protein
tmpProtein = Protein { 
    variance = tmp_variance,
    protein  = tmp_protein,
    lambda   = tmp_lambda
    }

-- | Пары, определяющие, в каком месте белка и на какую аминокислоту
-- мы можем произвести замену
bros_variance  = fst4 $ unzip4 bros_list
bros_position  = snd4 $ unzip4 bros_list
bros_prob_cros = trd4 $ unzip4 bros_list
bros_prob_mut  = fth4 $ unzip4 bros_list

-- | Вставка изменяемых аминокислот в шаблонный белок
-- с целью получить вид полученного белка  
insertVariance :: [(Aminoacid, Int)] -> [Aminoacid]
insertVariance x = foldl insert (protein tmpProtein) x
    where insert p (a, n) = take (n-1) p <> [a] <> drop n p

-- | Взятие случайной аминокислоты из набора доступных
selectAminoacid :: [Aminoacid] -> IO Aminoacid
selectAminoacid x = do
    r <- randomRIO (0, length x - 1)
    return $ x !! r