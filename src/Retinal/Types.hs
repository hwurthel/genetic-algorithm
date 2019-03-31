{-# LANGUAGE TemplateHaskell, TypeOperators, RankNTypes #-}
module Retinal.Types where


import Control.Category
import Data.Label
import Prelude hiding ((.), id)
import Text.Read (readMaybe)
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
newMolecule = mempty

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
        (x', y', z') = (get (x . coordin) atom, get (y . coordin) atom, get (z . coordin) atom)
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
-- Вовзращает возможные варианты молекул со вставкой
setFirstAtom :: ZMatrix -> Molecule -> [Molecule]
setFirstAtom zmatrix molecule =
     let matrixFirstAtom = zmatrix !! 1
         firstAtomID = fromJust $ get atomid   matrixFirstAtom
         zeroAtomID  = fromJust $ get atomcon  matrixFirstAtom
         distance    = fromJust $ get bonddist matrixFirstAtom

         firstAtom = get atom matrixFirstAtom
         zeroAtom' = molecule Map.!? zeroAtomID
         zeroAtom  = if isNothing zeroAtom' then error "setFirstAtom: zeroAtom not found" else fromJust zeroAtom'

         (fx, fy, fz) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) zeroAtom
         r1' = Point (fx - 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
         r2' = Point (fx + 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
         r3' = Point (fx - 2 * distance) (fy + 2 * distance) (fz - 2 * distance)
         r4' = Point (fx - 2 * distance) (fy - 2 * distance) (fz + 2 * distance)
         ws = sortCoordinates [r1',r2',r3',r4']

         alpha = [0, 1 .. 360]; beta  = [0, 1 .. 180]; 
         radians t = t * pi / 180

         -- | Вставляем @firstAtom@ в @molecule@
         -- | Отбираем только те молекулы, который находятся в окрестности @firstAtom@.
         -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
         -- с которым связан прикрепляемый атом (то есть удаляем @firstAtom@)
         wsMolecule = Map.elems $ Map.delete zeroAtomID $ Map.filter (`isInWorkSpace` ws) molecule

         -- | Функция, задающая предполагаемые координаты атома @tAtom@
         -- по углам радиусу, азимутальному и полярному углам (@dist@, @a@, @b@)
         updateCoord atom a b d = 
            let x' = fx + d*cos(radians a)*sin(radians b)
                y' = fy + d*sin(radians a)*sin(radians b)
                z' = fz + d*cos(radians b)
            in  set (z . coordin) z' $ set (y . coordin) y' $ set (x . coordin) x' atom

         -- | Всевозможные варианты расположений атома @firstAtom@ вокруг @zeroAtom@
         allVariance  = (updateCoord firstAtom) <$> alpha <*> beta <*> pure distance
         goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
     in  map (\x -> Map.insert firstAtomID x molecule) goodVariance


-- | Функция вставляет второй атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой
setSecondtAtom :: ZMatrix -> Molecule -> [Molecule]
setSecondtAtom zmatrix molecule =
          let matrixFirstAtom = zmatrix !! 1
              firstAtomID = fromJust $ get atomid   matrixFirstAtom
              zeroAtomID  = fromJust $ get atomcon  matrixFirstAtom
              distance    = fromJust $ get bonddist matrixFirstAtom
     
              firstAtom = get atom matrixFirstAtom
              zeroAtom' = molecule Map.!? zeroAtomID
              zeroAtom  = if isNothing zeroAtom' then error "setFirstAtom: zeroAtom not found" else fromJust zeroAtom'
     
              (fx, fy, fz) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) zeroAtom
              r1' = Point (fx - 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
              r2' = Point (fx + 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
              r3' = Point (fx - 2 * distance) (fy + 2 * distance) (fz - 2 * distance)
              r4' = Point (fx - 2 * distance) (fy - 2 * distance) (fz + 2 * distance)
              ws = sortCoordinates [r1',r2',r3',r4']
     
              alpha = [0, 1 .. 360]; beta  = [0, 1 .. 180]; 
              radians t = t * pi / 180
     
              -- | Вставляем @firstAtom@ в @molecule@
              -- | Отбираем только те молекулы, который находятся в окрестности @firstAtom@.
              -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
              -- с которым связан прикрепляемый атом (то есть удаляем @firstAtom@)
              wsMolecule = Map.elems $ Map.delete zeroAtomID $ Map.filter (`isInWorkSpace` ws) molecule
     
              -- | Функция, задающая предполагаемые координаты атома @tAtom@
              -- по углам радиусу, азимутальному и полярному углам (@dist@, @a@, @b@)
              updateCoord atom a b d = 
                 let x' = fx + d*cos(radians a)*sin(radians b)
                     y' = fy + d*sin(radians a)*sin(radians b)
                     z' = fz + d*cos(radians b)
                 in  set (z . coordin) z' $ set (y . coordin) y' $ set (x . coordin) x' atom
     
              -- | Всевозможные варианты расположений атома @firstAtom@ вокруг @zeroAtom@
              allVariance  = (updateCoord firstAtom) <$> alpha <*> beta <*> pure distance
              goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
          in  map (\x -> Map.insert firstAtomID x molecule) goodVariance

          setFirstAtom :: ZMatrix -> Molecule -> [Molecule]
          setFirstAtom zmatrix molecule =
               let matrixFirstAtom = zmatrix !! 1
                   firstAtomID = fromJust $ get atomid   matrixFirstAtom
                   zeroAtomID  = fromJust $ get atomcon  matrixFirstAtom
                   distance    = fromJust $ get bonddist matrixFirstAtom
          
                   firstAtom = get atom matrixFirstAtom
                   zeroAtom' = molecule Map.!? zeroAtomID
                   zeroAtom  = if isNothing zeroAtom' then error "setFirstAtom: zeroAtom not found" else fromJust zeroAtom'
          
                   (fx, fy, fz) = (\a -> (get (x . coordin) a, get (y . coordin) a, get (z . coordin) a)) zeroAtom
                   r1' = Point (fx - 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
                   r2' = Point (fx + 2 * distance) (fy - 2 * distance) (fz - 2 * distance)
                   r3' = Point (fx - 2 * distance) (fy + 2 * distance) (fz - 2 * distance)
                   r4' = Point (fx - 2 * distance) (fy - 2 * distance) (fz + 2 * distance)
                   ws = sortCoordinates [r1',r2',r3',r4']
          
                   alpha = [0, 1 .. 360]; beta  = [0, 1 .. 180]; 
                   radians t = t * pi / 180
          
                   -- | Вставляем @firstAtom@ в @molecule@
                   -- | Отбираем только те молекулы, который находятся в окрестности @firstAtom@.
                   -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
                   -- с которым связан прикрепляемый атом (то есть удаляем @firstAtom@)
                   wsMolecule = Map.elems $ Map.delete zeroAtomID $ Map.filter (`isInWorkSpace` ws) molecule
          
                   -- | Функция, задающая предполагаемые координаты атома @tAtom@
                   -- по углам радиусу, азимутальному и полярному углам (@dist@, @a@, @b@)
                   updateCoord atom a b d = 
                      let x' = fx + d*cos(radians a)*sin(radians b)
                          y' = fy + d*sin(radians a)*sin(radians b)
                          z' = fz + d*cos(radians b)
                      in  set (z . coordin) z' $ set (y . coordin) y' $ set (x . coordin) x' atom
          
                   -- | Всевозможные варианты расположений атома @firstAtom@ вокруг @zeroAtom@
                   allVariance  = (updateCoord firstAtom) <$> alpha <*> beta <*> pure distance
                   goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
               in  map (\x -> Map.insert firstAtomID x molecule) goodVariance
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
readMolecule inp = foldr addAtom Map.empty atoms
    where s = lines . unsafePerformIO $ readFile inp
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