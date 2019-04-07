{-# LANGUAGE TemplateHaskell  #-}
{-# LANGUAGE BangPatterns     #-}

module InsertMolecule
    ( setAtomWithOutOptimization
    , setAtomWithOptimization
    , readMolecule
    , readZMatrix
    , writeMolecule
    ) where

import Control.DeepSeq
import Control.Category
import Data.Label
import Prelude hiding ((.), id)
import Text.Read (readMaybe)
import Text.Printf
import Data.Maybe
import System.IO
import System.Directory
import System.IO.Unsafe
import Control.Concurrent

import qualified Data.Map as Map
import qualified Data.List as List
import qualified System.IO.Strict as StrictIO

import Utils

-- | НАЧАЛО. ОПИСАНИЕ ТИПОВ.
type CompType = Float
type ID      = Int
type Point   = (CompType, CompType, CompType)
type Name    = String
type Element = String
type Radius  = CompType

data ZElement = ZElement { _atom      :: Atom
                         , _atomid    :: Maybe ID
                         , _atomcon   :: Maybe ID
                         , _bonddist  :: Maybe CompType
                         , _anglcon   :: Maybe ID
                         , _bondangl  :: Maybe CompType
                         , _dihedcon  :: Maybe ID
                         , _dihedangl :: Maybe CompType
                         } deriving Show
type ZMatrix  = [ZElement]

data Atom = Atom { _name     :: Name
                 , _resname  :: Name
                 , _resseq    :: ID
                 , _coordin  :: Point
                 , _element  :: Element
                 , _radius   :: Radius
                 } deriving Show
type Molecule = Map.Map ID Atom

mkLabels [''Atom, ''ZElement]
-- | КОНЕЦ. ОПИСАНИЕ ТИПОВ.

-- | НАЧАЛО. КОНСТРУКТОРЫ.
-- Простейшие конструкторы пустых значений типов
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
newAtom = Atom { _name    = []
               , _resname = []
               , _resseq  = 0
               , _coordin = (0, 0, 0)
               , _element = []
               , _radius  = 0
               }

newMolecule :: Molecule
newMolecule = Map.empty

addAtom :: (ID, Atom) -> Molecule -> Molecule
addAtom (id, atom) molecule = Map.insert id atom molecule
-- | КОНЕЦ. КОНСТРУКТОРЫ.

-- | НАЧАЛО. ВСТАВКА РЕТИНАЛЯ.
-- | Функция предназначена для последовательной вставки
-- @n@ атомов молекулы из @zmatr@, находящихся в строках [1, n],
-- без процесса оптимизации.
-- В качестве результата возвращается первая молекулу, в которую удалось
-- вставить @n@ атомов. Если такой молекулы не нашлось, то возвращается Nothing
setAtomWithOutOptimization :: Int -> ZMatrix -> Molecule -> Maybe Molecule
setAtomWithOutOptimization n zmatr mol =
    listToMaybe $ setAtomWithOutOptimization' 1 n zmatr mol
    where setAtomWithOutOptimization' m n zmatr mol
            | m > n     = [mol]
            | m == 1    = concat $! setAtomWithOutOptimization' (m+1) n zmatr <$> setAtomWithOutOptimization1 1 zmatr mol
            | m == 2    = concat $! setAtomWithOutOptimization' (m+1) n zmatr <$> setAtomWithOutOptimization2 2 zmatr mol
            | otherwise = concat $! setAtomWithOutOptimization' (m+1) n zmatr <$> setAtomWithOutOptimization3 m zmatr mol

-- | Функция предназначена для последовательной вставки
-- атомов молекулы из @zmatr@, находящихся в строках [s, e],
-- c процессом оптимизации. В качестве результата возвращается молекула после всех вставок
setAtomWithOptimization :: Int -> Int -> ZMatrix -> Molecule -> IO Molecule
setAtomWithOptimization s e zmatr mol
    | s > e = return mol
    | otherwise = do mol' <- setAtomWithOptimization3 s zmatr mol
                     mol' `seq` setAtomWithOptimization (s+1) e zmatr mol'

                     -- | Функция предназначена для вставки
-- первого атома молекулы из @zmatrix@ без процесса оптимизации.
-- Возвращает ВСЕ возможные варианты молекул со вставкой.
-- Если таковых нет, то возвращается пустой список.
setAtomWithOutOptimization1 :: Int -> ZMatrix -> Molecule -> [Molecule]
setAtomWithOutOptimization1 n zmatrix molecule =
    do
    let matrixAtom_B = zmatrix !! n
        atomID_B     = fromJust $ get atomid   matrixAtom_B
        atomID_A     = fromJust $ get atomcon  matrixAtom_B
        distance_AB  = fromJust $ get bonddist matrixAtom_B

        atom_B  = get atom matrixAtom_B
        atom_A  = fromMaybe (error "setFirstAtom: atom_A not found") $ molecule Map.!? atomID_A

        (x_A, y_A, z_A) = get coordin atom_A

        -- | Отбираем только те молекулы, который находятся в окрестности @atom_A@.
        -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
        -- с которым связан прикрепляемый атом (то есть удаляем @atom_A@)
        m  = 2
        r1 = (,,) (x_A - m * distance_AB) (y_A - m * distance_AB) (z_A - m * distance_AB)
        r2 = (,,) (x_A + m * distance_AB) (y_A - m * distance_AB) (z_A - m * distance_AB)
        r3 = (,,) (x_A - m * distance_AB) (y_A + m * distance_AB) (z_A - m * distance_AB)
        r4 = (,,) (x_A - m * distance_AB) (y_A - m * distance_AB) (z_A + m * distance_AB)
        ws = sortCoordinates [r1,r2,r3,r4]
        wsMolecule = Map.elems $ Map.filter (`isInWorkSpace` ws) $ Map.delete atomID_A molecule

        -- | Вставляем @atom_B@ в @molecule@
        -- | Функция, задающая координаты атома @atom_B@ в единой СК
        -- по углам alpha и beta в штрихованной СК.
        possibleCoord a b =
            let x_B = distance_AB*cosd(a)*sind(b) + x_A
                y_B = distance_AB*sind(a)*sind(b) + y_A
                z_B = distance_AB*cosd(b)         + z_A
            in  (x_B, y_B, z_B) `deepseq` set coordin (x_B, y_B, z_B) atom_B
    let alpha = [Degree 0, Degree 5 .. Degree 375]
        beta  = [Degree 0, Degree 5 .. Degree 175]
        allVariance = possibleCoord <$> alpha <*> beta
    goodVariance <- filter (\x -> not . or $ isIntersection <$> wsMolecule <*> pure x) allVariance
    return $ Map.insert atomID_B goodVariance molecule

-- | Функция предназначена для вставки
-- второго атома молекулы из @zmatrix@ без процесса оптимизации.
-- Вовзращает ВСЕ возможные варианты молекул со вставкой.
-- Если таковых нет, то возвращается пустой список.
setAtomWithOutOptimization2 :: Int -> ZMatrix -> Molecule -> [Molecule]
setAtomWithOutOptimization2 n zmatrix molecule =
    do
    let matrixAtom_C = zmatrix !! n
        atomID_C     = fromJust $ get atomid   matrixAtom_C
        atomID_B     = fromJust $ get atomcon  matrixAtom_C
        atomID_A     = fromJust $ get anglcon  matrixAtom_C
        distance_BC  = fromJust $ get bonddist matrixAtom_C
        angle_ABC    = toDegree $ fromJust $ get bondangl matrixAtom_C

        atom_C  = get atom matrixAtom_C
        atom_B  = fromMaybe (error "setSecondAtom: atom_B not found") $ molecule Map.!? atomID_B
        atom_A  = fromMaybe (error "setSecondAtom: atom_A not found") $ molecule Map.!? atomID_A

        (x_A, y_A, z_A) = get coordin atom_A
        (x_B, y_B, z_B) = get coordin atom_B
        distance_AB = sqrt $ (x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2

        -- | Отбираем только те молекулы, которые находятся в окрестности @atom_B@.
        -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
        -- с которым связан прикрепляемый атом (то есть удаляем @atom_B@)
        m  = 2
        r1 = (,,) (x_B - m * distance_BC) (y_B - m * distance_BC) (z_B - m * distance_BC)
        r2 = (,,) (x_B + m * distance_BC) (y_B - m * distance_BC) (z_B - m * distance_BC)
        r3 = (,,) (x_B - m * distance_BC) (y_B + m * distance_BC) (z_B - m * distance_BC)
        r4 = (,,) (x_B - m * distance_BC) (y_B - m * distance_BC) (z_B + m * distance_BC)
        ws = sortCoordinates [r1,r2,r3,r4]
        wsMolecule = Map.elems $ Map.filter (`isInWorkSpace` ws) $ Map.delete atomID_B molecule

        -- | Ищем углы трансляции. Система координат левая.
        -- Поэтому положительным углам соответствует вращение по часовой стрелке.
        -- b1 - угол поворота вокруг оси Z
        -- b2 - угол поворота вокруг оси Y
        b1
            | xy == 0   =  Degree 0
            | y'_B > 0  =  acosd (x'_B / xy)
            | otherwise = -acosd (x'_B / xy)
            where
                xy = sqrt $ x'_B^2 + y'_B^2
                (x'_B, y'_B) = (x_B - x_A, y_B - y_A)

        b2
            | z'_B < 0  =  acosd (xy / xyz)
            | otherwise = -acosd (xy / xyz)
            where
                xy  = sqrt $ x'_B^2 + y'_B^2
                xyz = distance_AB
                (x'_B, y'_B, z'_B) = (x_B - x_A, y_B - y_A, z_B - z_A)

        -- | Вставляем @atom_C@ в @molecule@
        -- | Функция, задающая координаты атома @atom_C@  в единой СК
        -- по углу alpha в дважды штрихованной СК.
        possibleCoord a =
            let h     =  distance_BC * sind (angle_ABC)
                x''_C = -distance_BC * cosd (angle_ABC) + distance_AB
                z''_C =  h * sind (a)
                y''_C =  h * cosd (a)
                x_C   =  x''_C*cosd(b1)*cosd(b2) - y''_C*sind(b1) + z''_C*cosd(b1)*sind(b2) + x_A
                y_C   =  x''_C*sind(b1)*cosd(b2) + y''_C*cosd(b1) + z''_C*sind(b1)*sind(b2) + y_A
                z_C   = -x''_C*sind(b2)          + 0              + z''_C*cosd(b2)          + z_A
            in  (x_C, y_C, z_C) `seq` set coordin (x_C, y_C, z_C) atom_C
    let alpha = [Degree 0, Degree 2 .. Degree 358]
        allVariance  = possibleCoord <$> alpha
    goodVariance <- filter (\x -> not . or $ isIntersection <$> wsMolecule <*> pure x) allVariance
    return $ Map.insert atomID_C goodVariance molecule

-- | Функция предназначена для вставки
-- третьего и всех последующих атомов молекулы из @zmatrix@ без процесса оптимизации.
-- Вовзращает ОДИН возможный вариант молекулы со вставкой.
-- Если такового нет, то возвращается пустой список.
setAtomWithOutOptimization3 :: Int -> ZMatrix -> Molecule -> [Molecule]
setAtomWithOutOptimization3 n zmatrix molecule =
    do
    let matrixAtom_D = zmatrix !! n
        atomID_D     = fromJust $ get atomid   matrixAtom_D
        atomID_C     = fromJust $ get atomcon  matrixAtom_D
        atomID_B     = fromJust $ get anglcon  matrixAtom_D
        atomID_A     = fromJust $ get dihedcon matrixAtom_D
        distance_CD  = fromJust $ get bonddist matrixAtom_D
        angle_BCD    = toDegree $ fromJust $ get bondangl  matrixAtom_D
        angle_ABCD   = toDegree $ fromJust $ get dihedangl matrixAtom_D

        atom_D  = get atom matrixAtom_D
        atom_C  = fromMaybe (error "setThirdAtom: atom_C not found") $ molecule Map.!? atomID_C
        atom_B  = fromMaybe (error "setThirdAtom: atom_B not found") $ molecule Map.!? atomID_B
        atom_A  = fromMaybe (error "setThirdAtom: atom_A not found") $ molecule Map.!? atomID_A

        (x_A, y_A, z_A) = get coordin atom_A
        (x_B, y_B, z_B) = get coordin atom_B
        (x_C, y_C, z_C) = get coordin atom_C
        distance_BC = sqrt $ (x_C - x_B)^2 + (y_C - y_B)^2 + (z_C - z_B)^2

        -- | Отбираем только те молекулы, который находятся в окрестности @atom_B@.
        -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
        -- с которым связан прикрепляемый атом (то есть удаляем @atom_B@)
        -- !!! На самом деле здесь следует удалять из рассмотреия все атомы, с которыми связан вставлямый атом.
        m  = 2
        r1 = (,,) (x_B - m * distance_CD) (y_B - m * distance_CD) (z_B - m * distance_CD)
        r2 = (,,) (x_B + m * distance_CD) (y_B - m * distance_CD) (z_B - m * distance_CD)
        r3 = (,,) (x_B - m * distance_CD) (y_B + m * distance_CD) (z_B - m * distance_CD)
        r4 = (,,) (x_B - m * distance_CD) (y_B - m * distance_CD) (z_B + m * distance_CD)
        ws = sortCoordinates [r1,r2,r3,r4]
        wsMolecule = Map.elems $ Map.filter (`isInWorkSpace` ws) $ Map.delete atomID_C molecule

        -- | Ищем углы трансляции. Система координат левая. Упорядоченная тройка (x,y,z).
        -- Поэтому положительным углам соответствует вращение по часовой стрелке.
        -- b1 - угол поворота вокруг оси Z
        -- b2 - угол поворота вокруг оси Y
        b1  | xy == 0  =   Degree 0
            | y'_C > 0 =   acosd (x'_C / xy)
            | otherwise = -acosd (x'_C / xy)
            where
                xy = sqrt $ x'_C^2 + y'_C^2
                (x'_C, y'_C) = (x_C - x_B, y_C - y_B)

        b2  | z'_C < 0 =   acosd (xy / xyz)
            | otherwise = -acosd (xy / xyz)
            where
                xy  = sqrt $ x'_C^2 + y'_C^2
                xyz = distance_BC
                (x'_C, y'_C, z'_C) = (x_C - x_B, y_C - y_B, z_C - z_B)

        -- | Найдем координаты @atom_A@ в дважды штрихованной
        -- системе координат
        x'_A  =  x_A - x_B
        y'_A  =  y_A - y_B
        z'_A  =  z_A - z_B
        y''_A = -x'_A*sind(b1)          + y'_A*cosd(b1)          + 0
        z''_A =  x'_A*cosd(b1)*sind(b2) + y'_A*sind(b1)*sind(b2) + z'_A*cosd(b2)

        -- | Определеим угол между осью Z'' и направлением
        -- связи AB (угол отсчитывается в направлении по часовой стрелке)
        b3  | yz == 0 =    Degree 0
            | y''_A < 0 =  acosd (z''_A / yz)
            | otherwise = -acosd (z''_A / yz)
            where yz = sqrt $ y''_A^2 + z''_A^2

        -- | Вставляем @atom_C@ в @molecule@
        -- | Функция, задающая координаты атома @atom_C@  в единой СК
        -- в дважды штрихованной СК.
        possibleCoord =
            let h     =  distance_BC * sind (angle_BCD)
                x''_D = -distance_CD * cosd (angle_BCD) + distance_BC
                y''_D = -h * sind (b3 + angle_ABCD)
                z''_D =  h * cosd (b3 + angle_ABCD)
                x_D   =  x''_D*cosd(b1)*cosd(b2) - y''_D*sind(b1) + z''_D*cosd(b1)*sind(b2) + x_B
                y_D   =  x''_D*sind(b1)*cosd(b2) + y''_D*cosd(b1) + z''_D*sind(b1)*sind(b2) + y_B
                z_D   = -x''_D*sind(b2)          + 0              + z''_D*cosd(b2)          + z_B
            in  (x_D, y_D, z_D) `seq` set coordin (x_D, y_D, z_D) atom_D
    let allVariance = pure possibleCoord
    goodVariance <- filter (\x -> not . or $ isIntersection <$> wsMolecule <*> pure x) allVariance
    return $ Map.insert atomID_D possibleCoord molecule

-- | Функция предназначена для вставки
-- третьего и всех последующих атомов молекулы из @zmatrix@ с процесом оптимизации.
-- Вовзращает ОДИН возможный вариант молекулы со вставкой.
setAtomWithOptimization3 :: Int -> ZMatrix -> Molecule -> IO Molecule
setAtomWithOptimization3 n zmatrix molecule =
    do
    let matrixAtom_D = zmatrix !! n
        atomID_D     = fromJust $ get atomid   matrixAtom_D
        atomID_C     = fromJust $ get atomcon  matrixAtom_D
        atomID_B     = fromJust $ get anglcon  matrixAtom_D
        atomID_A     = fromJust $ get dihedcon matrixAtom_D
        distance_CD  = fromJust $ get bonddist matrixAtom_D
        angle_BCD    = toDegree $ fromJust $ get bondangl  matrixAtom_D
        angle_ABCD   = toDegree $ fromJust $ get dihedangl matrixAtom_D

        atom_D  = get atom matrixAtom_D
        atom_C  = fromMaybe (error "setAtomWithOptimization3: atom_C not found") $ molecule Map.!? atomID_C
        atom_B  = fromMaybe (error "setAtomWithOptimization3: atom_B not found") $ molecule Map.!? atomID_B
        atom_A  = fromMaybe (error "setAtomWithOptimization3: atom_A not found") $ molecule Map.!? atomID_A

        (x_A, y_A, z_A) = get coordin atom_A
        (x_B, y_B, z_B) = get coordin atom_B
        (x_C, y_C, z_C) = get coordin atom_C
        distance_BC = sqrt $ (x_C - x_B)^2 + (y_C - y_B)^2 + (z_C - z_B)^2

        -- | Ищем углы трансляции. Система координат левая. Упорядоченная тройка (x,y,z).
        -- Поэтому положительным углам соответствует вращение по часовой стрелке.
        -- b1 - угол поворота вокруг оси Z
        -- b2 - угол поворота вокруг оси Y
        b1  | xy == 0  =   Degree 0
            | y'_C > 0 =   acosd (x'_C / xy)
            | otherwise = -acosd (x'_C / xy)
                where
                    xy = sqrt $ x'_C^2 + y'_C^2
                    (x'_C, y'_C) = (x_C - x_B, y_C - y_B)

        b2  | z'_C < 0 =   acosd (xy / xyz)
            | otherwise = -acosd (xy / xyz)
            where
                xy  = sqrt $ x'_C^2 + y'_C^2
                xyz = distance_BC
                (x'_C, y'_C, z'_C) = (x_C - x_B, y_C - y_B, z_C - z_B)

        -- | Найдем координаты @atom_A@ в дважды штрихованной
        -- системе координат
        x'_A  =  x_A - x_B
        y'_A  =  y_A - y_B
        z'_A  =  z_A - z_B
        y''_A = -x'_A*sind(b1)          + y'_A*cosd(b1)          + 0
        z''_A =  x'_A*cosd(b1)*sind(b2) + y'_A*sind(b1)*sind(b2) + z'_A*cosd(b2)

        -- | Определеим угол между осью Z и направлением
        -- связи AB (угол отсчитывается в направлении по часовой стрелке)
        b3  | yz == 0   =  Degree 0
            | y''_A < 0 =  acosd (z''_A / yz)
            | otherwise = -acosd (z''_A / yz)
            where yz = sqrt $ y''_A^2 + z''_A^2

        -- | Вставляем @atom_C@ в @molecule@
        -- | Функция, задающая координаты атома @atom_C@  в единой СК
        -- по углу alpha в дважды штрихованной СК.
        possibleCoord =
            let h     =  distance_BC * sind (angle_BCD)
                x''_D = -distance_CD * cosd (angle_BCD) + distance_BC
                y''_D = -h * sind (b3 + angle_ABCD)
                z''_D =  h * cosd (b3 + angle_ABCD)
                x_D   =  x''_D*cosd(b1)*cosd(b2) - y''_D*sind(b1) + z''_D*cosd(b1)*sind(b2) + x_B
                y_D   =  x''_D*sind(b1)*cosd(b2) + y''_D*cosd(b1) + z''_D*sind(b1)*sind(b2) + y_B
                z_D   = -x''_D*sind(b2)          + 0              + z''_D*cosd(b2)          + z_B
            in  set coordin (x_D, y_D, z_D) atom_D
    let molecule' = Map.insert atomID_D possibleCoord molecule
        (x_D, y_D, z_D) = get coordin possibleCoord
        (r_D, m) = (get radius possibleCoord, 1)
        r1 = (,,) (x_D - m * r_D) (y_D - m * r_D) (z_D - m * r_D)
        r2 = (,,) (x_D + m * r_D) (y_D - m * r_D) (z_D - m * r_D)
        r3 = (,,) (x_D - m * r_D) (y_D + m * r_D) (z_D - m * r_D)
        r4 = (,,) (x_D - m * r_D) (y_D - m * r_D) (z_D + m * r_D)
        ws = sortCoordinates [r1,r2,r3,r4]
        wsMolecule     = Map.elems $ Map.filter (`isInWorkSpace` ws) molecule'
        atomsForOptim  = filter (isIntersection possibleCoord) wsMolecule
        resSeqForOptim = List.delete (get resseq possibleCoord) . List.nub $ map (get resseq) atomsForOptim
    let ouf = "for_optimization"
        inf = "optimizated_molecule"
        time_wait = 5000000
        wait False = return ()
        wait True  = do 
            e <- doesFileExist inf 
            if e then wait False
            else threadDelay time_wait >> wait True
    writeFile ouf $ unlines $ show <$> resSeqForOptim
    wait True
    molecule <- readMolecule inf
    removeFile ouf
    removeFile inf
    return molecule

-- | Функция сортирует координаты так, чтобы вектора
-- r0r1, r0r2, r0r3 образовывали правую тройку.
sortCoordinates :: [Point] -> [Point]
sortCoordinates points =
    let compX (x1,_,_) (x2,_,_) = compare x1 x2
        compY (_,y1,_) (_,y2,_) = compare y1 y2
        compZ (_,_,z1) (_,_,z2) = compare z1 z2
    in List.sortBy compZ . List.sortBy compY . List.sortBy compX $ points

-- | Функция определяет, находится ли атом в выделенном @[point]@ объеме.
isInWorkSpace :: Atom -> [Point] -> Bool
isInWorkSpace atom points =
    let r0:r1:r2:r3:_ = points
        (x0, y0, z0)  = r0
        (x1,  _,  _)  = r1
        ( _, y2,  _)  = r2
        ( _,  _, z3)  = r3
        (x', y', z') = get coordin atom
        between x a b = a <= x && x <= b
    in  and [between x' x0 x1, between y' y0 y2, between z' z0 z3]

-- | Функция принимает два атома и определяет, пересекаются
-- ли их Ван-дер-Вальсовы радиусы. В случае пересечения возвращается True.
isIntersection :: Atom -> Atom -> Bool
isIntersection atomA atomB =
    let rA = get radius atomA; (xA, yA, zA) = get coordin atomA
        rB = get radius atomB; (xB, yB, zB) = get coordin atomB
        dist = sqrt $ (xB - xA)^2 + (yB - yA)^2 + (zB - zA)^2
    in  dist <= (rA + rB)

-- | КОНЕЦ. ВСТАВКА РЕТИНАЛЯ.

-- | НАЧАЛО. ЧТЕНИЕ, ВЫВОД.
-- | Функция считывает молекулу из файла
readMolecule :: FilePath -> IO Molecule
readMolecule inf = return $ foldr addAtom newMolecule atoms
    where
        txt   = unsafePerformIO $ StrictIO.readFile inf
        atoms = [(id, Atom name resname resseq coordin elem radius) |
                 line <- lines txt,
                 let fields  = words line
                     id      = read $ fields !! 1
                     name    = fields !! 2
                     elem    = take 1 name
                     resname = fields !! 3
                     resseq   = read $ fields !! 5
                     coordin = (read $ fields !! 6, read $ fields !! 7, read $ fields !! 8)
                     radius  = getRadius elem]

writeMolecule :: FilePath -> Molecule -> IO ()
writeMolecule ouf (!molecule) =
    do  hdl <- openFile ouf WriteMode
        mapM_ (writeAtom hdl) (Map.toList molecule)
        hClose hdl
    where writeAtom hdl (id, atom) =
            let record  = "ATOM"
                serial  = id
                aname   = get name atom
                altLoc  = ' '
                altname = get resname atom
                chainid = 'A'
                resSeq  = get resseq atom
                icode   = ' '
                (x,y,z) = get coordin atom
                occupan = 0.0 :: Double
                tempfac = 0.0 :: Double
                pattern = "%-6s%5d %4s%1c%3s %1c%4d%1c   %8.3f%8.3f%8.3f%6.2f%6.2f\n"
            in hPrintf hdl pattern record serial aname altLoc altname chainid resSeq icode x y z occupan tempfac

-- | Функция считывает Z-матрицу из файла
readZMatrix :: FilePath -> IO ZMatrix
readZMatrix inf = return $ zmatrix
    where
        txt     = unsafePerformIO $ StrictIO.readFile inf
        zmatrix = [ZElement atom atomid atomcon bonddist anglcon bondangl dihedcon dihedangl |
                   line <- lines txt,
                   let fields = words line
                       atomid    = readMaybe $ fields !! 0
                       atomcon   = readMaybe $ fields !! 1
                       bonddist  = readMaybe $ fields !! 2
                       anglcon   = readMaybe $ fields !! 3
                       bondangl  = readMaybe $ fields !! 4
                       dihedcon  = readMaybe $ fields !! 5
                       dihedangl = readMaybe $ fields !! 6
                       name      = fields !! 7
                       element   = take 1 name
                       coordin   = (0, 0, 0)
                       radius    = getRadius element
                       resname   = fields !! 8
                       resid     = read $ fields !! 9
                       atom      = Atom name resname resid coordin element radius]

getRadius :: Element -> CompType
getRadius elem =
     case elem of "C" -> 1.100 -- // Заменено с 1.782
                  "H" -> 0.001 -- // Заменено с 0.200
                  "N" -> 1.100 -- // Заменено с 1.648
                  "O" -> 1.100 -- // Заменено с 1.510
                  "S" -> 1.100 -- // Заменено с 1.782
-- | КОНЕЦ. ЧТЕНИЕ, ВЫВОД.