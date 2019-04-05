{-# LANGUAGE TemplateHaskell #-}
module InsertMolecule  
    ( setFirstAtom
    , setSecondAtom
    , setThirdAtom
    , setOtherAtom
    , readMolecule
    , readZMatrix
    ) where

import Control.Category
import Data.Label
import Prelude hiding ((.), id)
import Text.Read (readMaybe)
import Text.Printf
import Data.Maybe
import Trigonometry

import qualified Data.Map as Map
import qualified Data.List as List
import System.IO.Unsafe

-- | НАЧАЛО. ОПИСАНИЕ ТИПОВ.
type ID      = Int
type Point   = (Double, Double, Double)
type Name    = String
type Element = String
type Radius  = Double

data ZElement = ZElement { _atom      :: Atom
                         , _atomid    :: Maybe ID
                         , _atomcon   :: Maybe ID
                         , _bonddist  :: Maybe Double
                         , _anglcon   :: Maybe ID
                         , _bondangl  :: Maybe Double
                         , _dihedcon  :: Maybe ID
                         , _dihedangl :: Maybe Double
                         } deriving Show
type ZMatrix  = [ZElement]

data Atom = Atom { _name     :: Name
                 , _coordin  :: Point
                 , _element  :: Element
                 , _radius   :: Radius
                 } deriving Show
type Molecule = Map.Map ID Atom

mkLabels [''Atom, ''ZElement]
-- | КОНЕЦ. ОПИСАНИЕ ТИПОВ.

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
               , _coordin = (0, 0, 0)
               , _element = []
               , _radius  = 0
               }

newMolecule :: Molecule
newMolecule = Map.empty

addAtom :: (ID, Atom) -> Molecule -> Molecule
addAtom (id, atom) molecule = Map.insert id atom molecule

-- | Функция сортирует координаты так, чтобы вектора
-- r0r1, r0r2, r0r3 образовывали правую тройку.
sortCoordinates :: [Point] -> [Point]
sortCoordinates points =
    let compX (x1,_,_) (x2,_,_) = compare x1 x2
        compY (_,y1,_) (_,y2,_) = compare y1 y2
        compZ (_,_,z1) (_,_,z2) = compare z1 z2
    in List.sortBy compZ . List.sortBy compY . List.sortBy compX $ points

-- | Функция определяет, находится ли атом в выделенном объеме.
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

-- | Функция принимает два атома и определяет, пересекаются не пересекаются
-- ли их Ван-дер-Вальсовы радиусы. Если нет, то возвращается True, иначе - False.
isIntersection :: Atom -> Atom -> Bool
isIntersection atomA atomB = 
    let rA = get radius atomA; (xA, yA, zA) = get coordin atomA
        rB = get radius atomB; (xB, yB, zB) = get coordin atomB
        dist = sqrt $ (xB - xA)^2 + (yB - yA)^2 + (zB - zA)^2
    in  dist <= (rA + rB)

-- | НАЧАЛО. ВСТАВКА РЕТИНАЛЯ.
-- | Функция вставляет первый атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой.
setFirstAtom :: ZMatrix -> Molecule -> [Molecule]
setFirstAtom zmatrix molecule =
     let matrixAtom_B = zmatrix !! 1
         atomID_B     = fromJust $ get atomid   matrixAtom_B
         atomID_A     = fromJust $ get atomcon  matrixAtom_B
         distance_AB  = fromJust $ get bonddist matrixAtom_B

         atom_B  = get atom matrixAtom_B
         atom_A' = molecule Map.!? atomID_A
         atom_A  = if isNothing atom_A' then error "setFirstAtom: atom_A not found" else fromJust atom_A'
         
         (x_A, y_A, z_A) = get coordin atom_A

         -- | Отбираем только те молекулы, который находятся в окрестности @atom_A@.
         -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
         -- с которым связан прикрепляемый атом (то есть удаляем @atom_A@)
         n  = 2
         r1 = (,,) (x_A - n * distance_AB) (y_A - n * distance_AB) (z_A - n * distance_AB)
         r2 = (,,) (x_A + n * distance_AB) (y_A - n * distance_AB) (z_A - n * distance_AB)
         r3 = (,,) (x_A - n * distance_AB) (y_A + n * distance_AB) (z_A - n * distance_AB)
         r4 = (,,) (x_A - n * distance_AB) (y_A - n * distance_AB) (z_A + n * distance_AB)
         ws = sortCoordinates [r1,r2,r3,r4]
         wsMolecule = Map.elems $ Map.delete atomID_A $ Map.filter (`isInWorkSpace` ws) molecule
         
         -- | Вставляем @atom_B@ в @molecule@
         -- | Функция, задающая координаты атома @atom_B@ в единой СК 
         -- по углам alpha и beta в штрихованной СК.
         possibleCoord a b = 
            let x_B = distance_AB*cosd(a)*sind(b) + x_A
                y_B = distance_AB*sind(a)*sind(b) + y_A
                z_B = distance_AB*cosd(b)         + z_A 
            in  set coordin (x_B, y_B, z_B) atom_B

         -- | Всевозможные варианты расположений атома @atom_B@ вокруг @atom_A@
         (alpha, beta) = ([Degree 0, Degree 1 .. Degree 359], [Degree 0, Degree 1 .. Degree 179])
         allVariance  = possibleCoord <$> alpha <*> beta
         goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
     in  map (\x -> Map.insert atomID_B x molecule) goodVariance

-- | Функция вставляет второй атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой.
setSecondAtom :: ZMatrix -> Molecule -> [Molecule]
setSecondAtom zmatrix molecule =
          let matrixAtom_C = zmatrix !! 2
              atomID_C     = fromJust $ get atomid   matrixAtom_C
              atomID_B     = fromJust $ get atomcon  matrixAtom_C
              atomID_A     = fromJust $ get anglcon  matrixAtom_C
              distance_BC  = fromJust $ get bonddist matrixAtom_C
              angle_ABC    = toDegree $ fromJust $ get bondangl matrixAtom_C

              atom_C  = get atom matrixAtom_C
              atom_B' = molecule Map.!? atomID_B
              atom_A' = molecule Map.!? atomID_A
              atom_B  = if isNothing atom_B' then error "setSecondAtom: atom_B not found" else fromJust atom_B'
              atom_A  = if isNothing atom_A' then error "setSecondAtom: atom_A not found" else fromJust atom_A'

              (x_A, y_A, z_A) = get coordin atom_A
              (x_B, y_B, z_B) = get coordin atom_B
              distance_AB = sqrt $ (x_B - x_A)^2 + (y_B - y_A)^2 + (z_B - z_A)^2

              -- | Отбираем только те молекулы, которые находятся в окрестности @atom_B@.
              -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
              -- с которым связан прикрепляемый атом (то есть удаляем @atom_B@)
              n  = 2
              r1 = (,,) (x_B - n * distance_BC) (y_B - n * distance_BC) (z_B - n * distance_BC)
              r2 = (,,) (x_B + n * distance_BC) (y_B - n * distance_BC) (z_B - n * distance_BC)
              r3 = (,,) (x_B - n * distance_BC) (y_B + n * distance_BC) (z_B - n * distance_BC)
              r4 = (,,) (x_B - n * distance_BC) (y_B - n * distance_BC) (z_B + n * distance_BC)
              ws = sortCoordinates [r1,r2,r3,r4]
              wsMolecule = Map.elems $ Map.delete atomID_B $ Map.filter (`isInWorkSpace` ws) molecule

              -- | Ищем углы трансляции. Система координат левая.
              -- Поэтому положительным углам соответствует вращение по часовой стрелке.
              -- b1 - угол поворота вокруг оси Z
              -- b2 - угол поворота вокруг оси Y 
              b1 = 
                if xy == 0 then Degree 0
                else if y'_B > 0 
                     then  acosd (x'_B / xy)
                     else -acosd (x'_B / xy)
                where
                    xy = sqrt $ x'_B^2 + y'_B^2
                    (x'_B, y'_B) = (x_B - x_A, y_B - y_A)
              
              b2 = 
                if z'_B < 0
                then  acosd (xy / xyz)
                else -acosd (xy / xyz)
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
                 in  set coordin (x_C, y_C, z_C) atom_C
     
              -- | Всевозможные варианты расположений атома @atom_C@
              alpha = [Degree 0, Degree 10 .. Degree 350]
              allVariance  = possibleCoord <$> alpha
              goodVariance = map fst $ filter (\(_,b) -> b == False) $ map (\x -> (x, or $ isIntersection <$> wsMolecule <*> pure x)) allVariance
          in  map (\x -> Map.insert atomID_C x molecule) goodVariance

-- | Функция вставляет третий атом ретиналя. 
-- Вовзращает возможные варианты молекул со вставкой.
--setThirdAtom :: ZMatrix -> Molecule -> [Molecule]
setThirdAtom zmatrix molecule =
    let matrixAtom_D = zmatrix !! 3
        atomID_D     = fromJust $ get atomid   matrixAtom_D
        atomID_C     = fromJust $ get atomcon  matrixAtom_D
        atomID_B     = fromJust $ get anglcon  matrixAtom_D
        atomID_A     = fromJust $ get dihedcon matrixAtom_D
        distance_CD  = fromJust $ get bonddist matrixAtom_D
        angle_BCD    = toDegree $ fromJust $ get bondangl  matrixAtom_D
        angle_ABCD   = toDegree $ fromJust $ get dihedangl matrixAtom_D

        atom_D  = get atom matrixAtom_D
        atom_C' = molecule Map.!? atomID_C
        atom_B' = molecule Map.!? atomID_B
        atom_A' = molecule Map.!? atomID_A
        atom_C  = if isNothing atom_C' then error "setThirdAtom: atom_C not found" else fromJust atom_C'
        atom_B  = if isNothing atom_B' then error "setThirdAtom: atom_B not found" else fromJust atom_B'
        atom_A  = if isNothing atom_A' then error "setThirdAtom: atom_A not found" else fromJust atom_A'

        (x_A, y_A, z_A) = get coordin atom_A
        (x_B, y_B, z_B) = get coordin atom_B
        (x_C, y_C, z_C) = get coordin atom_C
        distance_BC = sqrt $ (x_C - x_B)^2 + (y_C - y_B)^2 + (z_C - z_B)^2

        -- | Отбираем только те молекулы, который находятся в окрестности @atom_B@.
        -- Для дальнейшего рассмотрения удаляем из рассмотрения тот атом,
        -- с которым связан прикрепляемый атом (то есть удаляем @atom_B@)
        -- !!! На самом деле здесь следует удалять из рассмотреия все атомы, с которыми связан вставлямый атом.
        n  = 2
        r1 = (,,) (x_B - n * distance_CD) (y_B - n * distance_CD) (z_B - n * distance_CD)
        r2 = (,,) (x_B + n * distance_CD) (y_B - n * distance_CD) (z_B - n * distance_CD)
        r3 = (,,) (x_B - n * distance_CD) (y_B + n * distance_CD) (z_B - n * distance_CD)
        r4 = (,,) (x_B - n * distance_CD) (y_B - n * distance_CD) (z_B + n * distance_CD)
        ws = sortCoordinates [r1,r2,r3,r4]
        wsMolecule = Map.elems $ Map.delete atomID_B $ Map.filter (`isInWorkSpace` ws) molecule

        -- | Ищем углы трансляции. Система координат левая. Упорядоченная тройка (x,y,z).
        -- Поэтому положительным углам соответствует вращение по часовой стрелке.
        -- b1 - угол поворота вокруг оси Z
        -- b2 - угол поворота вокруг оси Y 
        b1 = 
            if xy == 0 then Degree 0
            else if y'_C > 0 
                then  acosd (x'_C / xy)
                else -acosd (x'_C / xy)
            where
                xy = sqrt $ x'_C^2 + y'_C^2
                (x'_C, y'_C) = (x_C - x_B, y_C - y_B)
        
        b2 = 
            if z'_C < 0
            then  acosd (xy / xyz)
            else -acosd (xy / xyz)
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
        b3 = 
            if yz == 0 then Degree 0
            else if y''_A < 0
            then  acosd (z''_A / yz)
            else -acosd (z''_A / yz)
            where 
            yz = sqrt $ y''_A^2 + z''_A^2

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
            in  set coordin (x_D, y_D, z_D) atom_D      
     in [] 
        --if isIntersection <$> wsMolecule <*> pure possibleCoord
    --    then possibleCoord
    --    else []
            --     then Map.empty
            --     else head $ map (\x -> Map.insert atomID_D x molecule) goodVariance

-- | Функция вставляет все оставшиеся атомы ретиналя. 
setOtherAtom :: ZMatrix -> Int -> Molecule -> IO Molecule
setOtherAtom zmatrix n molecule =
            let matrixAtom_D = zmatrix !! n
                atomID_D     = fromJust $ get atomid   matrixAtom_D
                atomID_C     = fromJust $ get atomcon  matrixAtom_D
                atomID_B     = fromJust $ get anglcon  matrixAtom_D
                atomID_A     = fromJust $ get dihedcon matrixAtom_D
                distance_CD  = fromJust $ get bonddist matrixAtom_D
                angle_BCD    = toDegree $ fromJust $ get bondangl  matrixAtom_D
                angle_ABCD   = toDegree $ fromJust $ get dihedangl matrixAtom_D
  
                atom_D  = get atom matrixAtom_D
                atom_C' = molecule Map.!? atomID_C
                atom_B' = molecule Map.!? atomID_B
                atom_A' = molecule Map.!? atomID_A
                atom_C  = if isNothing atom_C' then error "setThirdAtom: atom_C not found" else fromJust atom_C'
                atom_B  = if isNothing atom_B' then error "setThirdAtom: atom_B not found" else fromJust atom_B'
                atom_A  = if isNothing atom_A' then error "setThirdAtom: atom_A not found" else fromJust atom_A'
  
                (x_A, y_A, z_A) = get coordin atom_A
                (x_B, y_B, z_B) = get coordin atom_B
                (x_C, y_C, z_C) = get coordin atom_C
                distance_BC = sqrt $ (x_C - x_B)^2 + (y_C - y_B)^2 + (z_C - z_B)^2
  
                -- | Ищем углы трансляции. Система координат левая. Упорядоченная тройка (x,y,z).
                -- Поэтому положительным углам соответствует вращение по часовой стрелке.
                -- b1 - угол поворота вокруг оси Z
                -- b2 - угол поворота вокруг оси Y 
                b1 = 
                  if xy == 0 then Degree 0
                  else if y'_C > 0 
                       then  acosd (x'_C / xy)
                       else -acosd (x'_C / xy)
                  where
                      xy = sqrt $ x'_C^2 + y'_C^2
                      (x'_C, y'_C) = (x_C - x_B, y_C - y_B)
                
                b2 = 
                  if z'_C < 0
                  then  acosd (xy / xyz)
                  else -acosd (xy / xyz)
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
                b3 = 
                  if yz == 0 then Degree 0
                  else if y''_A < 0
                    then  acosd (z''_A / yz)
                    else -acosd (z''_A / yz)
                  where 
                    yz = sqrt $ y''_A^2 + z''_A^2

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
            in do
                let molecule' = Map.insert atomID_D possibleCoord molecule
                print "Optimization" -- Здесь мы выподим молекулу в файл для оптимизатора 
                return $ readMolecule "optimizated molecule" -- Читаем оптимизированную молекулу 

-- ТЕСТ
--   let mol = readMolecule "protein"
--   let matr = readZMatrix  "z-matrix"
--   let molecule = (setFirstAtom matr mol) !! 2
--   let res = setSecondAtom matr molecule
--   let f a = printf "ATOM   3370   CA REL U 212      %.3f  %.3f    %.3f  1.00  0.00      U    C\n" (get (x. coordin) a) (get (y.coordin) a) (get (z.coordin) a)
--   mapM_ f res

-- | КОНЕЦ. ВСТАВКА РЕТИНАЛЯ.

-- | НАЧАЛО. ЧТЕНИЕ И ОБРАБОТКА ФАЙЛОВ.
-- | Функция считывает молекулу из файла
readMolecule :: FilePath -> Molecule
readMolecule inf = foldr addAtom newMolecule atoms
    where
        txt   = unsafePerformIO $ readFile inf
        atoms = [(id, Atom name coordin elem radius) |  
                 line <- lines txt,
                 let fields  = words line 
                     id      = read   $ fields !! 1
                     name    = fields !! 2
                     elem    = take 1 name
                     coordin = (read $ fields !! 6, read $ fields !! 7, read $ fields !! 8)
                     radius  = getRadius elem]

-- | Функция считывает Z-матрицу из файла
readZMatrix :: FilePath -> ZMatrix
readZMatrix inf = zmatrix
    where
        txt     = unsafePerformIO $ readFile inf 
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
                       atom      = Atom name coordin element radius]

getRadius :: Element -> Double
getRadius elem = 
     case elem of "C" -> 1.282 -- // Заменено с 1.782
                  "H" -> 0.200 
                  "N" -> 1.248 -- // Заменено с 1.648
                  "O" -> 1.515
                  "S" -> 1.782
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