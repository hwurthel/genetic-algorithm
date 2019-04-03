{-# LANGUAGE TemplateHaskell, TypeOperators, RankNTypes, GeneralizedNewtypeDeriving #-}
module Retinal.Types where


import Control.Category
import Data.Label
import Prelude hiding ((.), id)
import Text.Read (readMaybe)
import Text.Printf
import Data.Maybe

import qualified Data.Map as Map
import qualified Data.List as List
import System.IO.Unsafe

type    ID      = Int
newtype Name    = Name    String deriving Show
newtype Element = Element String deriving Show
newtype Radius  = Radius  Double deriving Show
data    Point   = Point { _x :: Double
                        , _y :: Double
                        , _z :: Double
                        } deriving (Eq, Ord)

instance Show Point where
     show (Point x y z) = unwords $ map show [x,y,z]                        

data ZElement = ZElement { _atom      :: Atom
                         , _atomid    :: Maybe ID
                         , _atomcon   :: Maybe ID
                         , _bonddist  :: Maybe Double
                         , _anglcon   :: Maybe ID
                         , _bondangl  :: Maybe Double
                         , _dihedcon  :: Maybe ID
                         , _dihedangl :: Maybe Double
                         } deriving Show

data Atom = Atom { _name     :: Name
                 , _coordin  :: Point
                 , _element  :: Element
                 , _radius   :: Radius
                 } deriving Show

type Molecule = Map.Map ID Atom
type ZMatrix  = [ZElement]

mkLabels [''Point, ''Atom, ''ZElement]

--instance Show Atom where
--    show = \x -> show $ _coordin x

-- | НАЧАЛО. ТРИГОНОМЕТРИЧЕСКИЕ ФУНКЦИИ, ГРАДУСЫ.
-- Работа с тригонометрическими функциями,
-- аргументам которых являются величины в градуса
newtype Degree = Degree Double deriving (Num, Enum, Show, Eq, Ord, Fractional)
toDegree :: Double -> Degree
toDegree = Degree

rad2Deg :: Double -> Degree
rad2Deg = \x -> Degree (x * 180 / pi)

deg2Rad :: Degree -> Double
deg2Rad = \(Degree x) -> x * pi / 180

cosd :: Degree -> Double
cosd = cos . deg2Rad 

sind :: Degree -> Double
sind = sin . deg2Rad

acosd :: Double -> Degree
acosd = rad2Deg . acos

asind :: Double -> Degree
asind = rad2Deg . asin
-- | КОНЕЦ. ТРИГОНОМЕТРИЧЕСКИЕ ФУНКЦИИ, ГРАДУСЫ.

newZElement :: ZElement
newZElement = ZElement { _atom      = newAtom
                       , _atomid    = Nothing
                       , _atomcon   = Nothing
                       , _bonddist  = Nothing
                       , _anglcon   = Nothing
                       , _bondangl  = Nothing
                       , _dihedcon  = Nothing 
                       , _dihedangl = Nothing
                       }

newZMatrix :: ZMatrix
newZMatrix = mempty

newAtom :: Atom
newAtom = Atom { _name    = Name []
               , _coordin = Point 0 0 0
               , _element = Element []
               , _radius  = Radius 0
               }

newMolecule :: Molecule
newMolecule = Map.empty

addAtom :: (ID, Atom) -> Molecule -> Molecule
addAtom (id, atom) molecule = Map.insert id atom molecule

-- | Функция сортирует координаты так, чтобы вектора
-- r0r1, r0r2, r0r3 образовывали правую тройку.
sortCoordinates :: [Point] -> [Point]
sortCoordinates points =
    let compX a@(Point x1 _ _) b@(Point x2 _ _) = compare x1 x2
        compY a@(Point _ y1 _) b@(Point _ y2 _) = compare y1 y2
        compZ a@(Point _ _ z1) b@(Point _ _ z2) = compare z1 z2
    in List.sortBy compZ . List.sortBy compY . List.sortBy compX $ points

-- | Функция определяет, находится ли атом в выделенном объеме.
isInWorkSpace :: Atom -> [Point] -> Bool
isInWorkSpace atom points =
    let r0:r1:r2:r3:_ = points
        (x', y', z') = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) atom
        between x a b = a <= x && x <= b
    in  and [between x' (get x r0) (get x r1), between y' (get y r0) (get y r2), between z' (get z r0) (get z r3)]

-- | Функция принимает два атома и определяет, пересекаются не пересекаются
-- ли их Ван-дер-Вальсовы радиусы. Если нет, то возвращается True, иначе - False.
isIntersection :: Atom -> Atom -> Bool
isIntersection atomA atomB = 
    let getRXYZ = \a -> (get radius a, get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)
        (Radius rA, xA, yA, zA) = getRXYZ atomA
        (Radius rB, xB, yB, zB) = getRXYZ atomB
        dist = sqrt $ (xB - xA)^2 + (yB - yA)^2 + (zB - zA)^2
    in  dist <= (rA + rB)


-- | НАЧАЛО. ВСТАВКА РЕТИНАЛЯ.
-- | Функция вставляет первый атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой.
--setFirstAtom :: ZMatrix -> Molecule -> [Molecule]
setFirstAtom zmatrix molecule =
     let matrixAtom_B = zmatrix !! 1
         atomID_B     = fromJust $ get atomid   matrixAtom_B
         atomID_A     = fromJust $ get atomcon  matrixAtom_B
         distance_AB  = fromJust $ get bonddist matrixAtom_B

         atom_B  = get atom matrixAtom_B
         atom_A' = molecule Map.!? atomID_A
         atom_A  = if isNothing atom_A' then error "setFirstAtom: atom_A not found" else fromJust atom_A'
         
         (x_A, y_A, z_A) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) atom_A

         -- | Отбираем только те молекулы, который находятся в окрестности @atom_A@.
         -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
         -- с которым связан прикрепляемый атом (то есть удаляем @atom_A@)
         n  = 2
         r1 = Point (x_A - n * distance_AB) (y_A - n * distance_AB) (z_A - n * distance_AB)
         r2 = Point (x_A + n * distance_AB) (y_A - n * distance_AB) (z_A - n * distance_AB)
         r3 = Point (x_A - n * distance_AB) (y_A + n * distance_AB) (z_A - n * distance_AB)
         r4 = Point (x_A - n * distance_AB) (y_A - n * distance_AB) (z_A + n * distance_AB)
         ws = sortCoordinates [r1,r2,r3,r4]
         wsMolecule = Map.elems $ Map.delete atomID_A $ Map.filter (`isInWorkSpace` ws) molecule
         
         -- | Вставляем @atom_B@ в @molecule@
         -- | Функция, задающая координаты атома @atom_B@ в единой СК 
         -- по углам alpha и beta в штрихованной СК.
         possibleCoord a b = 
            let x_B = distance_AB*cosd(a)*sind(b) + x_A
                y_B = distance_AB*sind(a)*sind(b) + y_A
                z_B = distance_AB*cosd(b)         + z_A 
            in  set (z . coordin) z_B . set (y . coordin) y_B . set (x . coordin) x_B $ atom_B

         -- | Всевозможные варианты расположений атома @atom_B@ вокруг @atom_A@
         (alpha, beta) = ([Degree 0, Degree 1 .. Degree 359], [Degree 0, Degree 1 .. Degree 179])
         allVariance  = possibleCoord <$> alpha <*> beta
         goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
     in  map (\x -> Map.insert atomID_B x molecule) goodVariance


-- | Функция вставляет второй атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой.
--setSecondtAtom :: ZMatrix -> Molecule -> [Molecule]
setSecondAtom zmatrix molecule =
          let matrixAtom_C = zmatrix !! 2
              atomID_C     = fromJust $ get atomid   matrixAtom_C
              atomID_B     = fromJust $ get atomcon  matrixAtom_C
              atomID_A     = fromJust $ get anglcon  matrixAtom_C
              distance_BC  = fromJust $ get bonddist matrixAtom_C
              angle_CBA    = toDegree $ fromJust $ get bondangl matrixAtom_C

              atom_C  = get atom matrixAtom_C
              atom_B' = molecule Map.!? atomID_B
              atom_A' = molecule Map.!? atomID_A
              atom_B  = if isNothing atom_B' then error "setSecondAtom: atom_B not found" else fromJust atom_B'
              atom_A  = if isNothing atom_A' then error "setSecondAtom: atom_A not found" else fromJust atom_A'

              (x_A, y_A, z_A) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) atom_A
              (x_B, y_B, z_B) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) atom_B  -- (16.910,  35.413,    5.030)
              distance_AB = sqrt $ (x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2
              distance_AC = sqrt $ distance_AB^2 + distance_BC^2 - 2*distance_AB*distance_BC*cosd(angle_CBA)
              angle_BAC   = acosd ((distance_AC^2 + distance_AB^2 - distance_BC^2 ) / (2*distance_AC*distance_AB))

              -- | Отбираем только те молекулы, который находятся в окрестности @atom_B@.
              -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
              -- с которым связан прикрепляемый атом (то есть удаляем @atom_B@)
              -- !!! На самом деле здесь следует удалять из рассмотреия все атомы, с которыми связан вставлямый атом.
              n  = 2
              r1 = Point (x_B - n * distance_BC) (y_B - n * distance_BC) (z_B - n * distance_BC)
              r2 = Point (x_B + n * distance_BC) (y_B - n * distance_BC) (z_B - n * distance_BC)
              r3 = Point (x_B - n * distance_BC) (y_B + n * distance_BC) (z_B - n * distance_BC)
              r4 = Point (x_B - n * distance_BC) (y_B - n * distance_BC) (z_B + n * distance_BC)
              ws = sortCoordinates [r1,r2,r3,r4]
              wsMolecule = Map.elems $ Map.delete atomID_B $ Map.filter (`isInWorkSpace` ws) molecule

              -- | Ищем углы трансляции.
              -- b1 - угол поворота вокруг оси Z
              -- b2 - угол поворота вокруг оси Y 
              b1 = 
                if xy == 0 then Degree 0
                else acosd (x'_B / xy)
                where
                    xy = sqrt $ x'_B^2 + y'_B^2
                    (x'_B, y'_B) = (x_B - x_A, y_B - y_A)
              
              b2 = acosd (xy / xyz)
                where
                    xy  = sqrt $ x'_B^2 + y'_B^2
                    xyz = distance_AB
                    (x'_B, y'_B, z'_B) = (x_B - x_A, y_B - y_A, z_B - z_A)
              
              -- | Вставляем @atom_C@ в @molecule@
              -- | Функция, задающая координаты атома @atom_C@  в единой СК 
              -- по углу alpha в дважды штрихованной СК.
              possibleCoord a = 
                 let h     =  distance_BC * sind (angle_CBA)
                     x''_C = -distance_BC * cosd (angle_CBA) + distance_AB
                     z''_C =  h * cosd (a)
                     y''_C =  h * sind (a)
                     x_C   =  x''_C*cosd(b1)*cosd(b2) - y''_C*sind(b1) + z''_C*cosd(b1)*sind(b2) + x_A
                     y_C   =  x''_C*sind(b1)*cosd(b2) + y''_C*cosd(b1) + z''_C*sind(b1)*sind(b2) + y_A
                     z_C   = -x''_C*sind(b2)          + 0              + z''_C*cosd(b2)          + z_A
                 in  set (z . coordin) z_C . set (y . coordin) y_C . set (x . coordin) x_C $ atom_C
     
              -- | Всевозможные варианты расположений атома @atom_C@
              alpha = [Degree 0, Degree 1 .. Degree 359]
              allVariance  = possibleCoord <$> alpha
              goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
          in  map (\x -> Map.insert atomID_C x molecule) goodVariance

-- ТЕСТ
--   let mol = readMolecule "protein"
--   let matr = readZMatrix  "z-matrix"
--   let molecule = (setFirstAtom matr mol) !! 2
--   let res = setSecondAtom matr molecule
--   let f a = printf "ATOM   3370   CA REL U 212      %.3f  %.3f    %.3f  1.00  0.00      U    C\n" (get (x. coordin) a) (get (y.coordin) a) (get (z.coordin) a)
--   mapM_ f res

-- | КОНЕЦ. ВСТАВКА РЕТИНАЛЯ.

          -- | НАЧАЛО. ЧТЕНИЕ И ОБРАБОТКА ФАЙЛОВ.
-- | Функция преобразует файл @inf@ в формат, с которым работает программа
reWriteFile :: FilePath -> FilePath -> IO ()
reWriteFile inf ouf =  writeFile ouf (unlines newtxt)
     where s = unsafePerformIO $ readFile inf
           newtxt = [ unwords $ a0:a1:a2:a3:a4:a5:a6:[] |
                    line <- lines s,
                    let (_:a0:a1:_:_:_:a2:a3:a4:_:_:_:a5:[]) = words line,
                    let a6 = case a5 of   "C" -> "1.282" -- // Заменено с 1.782
                                          "H" -> "0.200"
                                          "N" -> "1.248" -- // Заменено с 1.648
                                          "O" -> "1.515"
                                          "S" -> "1.782"]

-- | Функция считывает молекулу из файла
readMolecule :: FilePath -> Molecule
readMolecule inp = foldr addAtom newMolecule atoms
    where 
        s = lines . unsafePerformIO $ readFile inp
        atoms = map readAtom s
        readAtom str = (a0', Atom a1' a2' a3' a4')
            where 
                a0:a1:a2:a3:a4:a5:a6:_ = words str
                a0' = read a0
                a1' = Name a1
                a2' = Point (read a2) (read a3) (read a4)
                a3' = Element a5
                a4' = Radius (read a6)

-- | Функция считывает Z-матрицу из файла
readZMatrix :: FilePath -> ZMatrix
readZMatrix inf = map readZElement (lines . unsafePerformIO $ readFile inf)
    where 
        readZElement str = do
            let a0:a1:a2:a3:a4:a5:a6:a7:a8:a9:_ = words str
                a0' = readMaybe a0
                a1' = readMaybe a1
                a2' = readMaybe a2
                a3' = readMaybe a3
                a4' = readMaybe a4
                a5' = readMaybe a5
                a6' = readMaybe a6
                a7' = a7
                a8' = a8
                a9' = read a9
                atom = Atom (Name a7') (Point 0 0 0) (Element a8') (Radius a9') 
            ZElement atom a0' a1' a2' a3' a4' a5' a6'
-- | КОНЕЦ. ЧТЕНИЕ И ОБРАБОТКА ФАЙЛОВ.
                 
                 
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