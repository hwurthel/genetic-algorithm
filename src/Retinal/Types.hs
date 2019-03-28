module Retinal.Types where

import qualified Data.Map as Map
import qualified Data.List as List
import System.IO.Unsafe

-- | Ортогональный координаты в ангстремах. 
data Point = Point 
    { getX :: Double
    , getY :: Double
    , getZ :: Double
    } deriving (Eq, Ord)

-- | Параметры атомов
newtype ID      = ID      Int    deriving (Show, Eq, Ord)
newtype Name    = Name    String deriving Show
newtype ResName = ResName String deriving Show
newtype ResSeq  = ResSeq  Int    deriving Show
newtype Element = Element String deriving Show
newtype Radius  = Radius  Double deriving Show

data Atom = Atom 
    { name     :: Name
    , resName  :: ResName
    , resSeq   :: ResSeq
    , coordin  :: Point
    , element  :: Element
    , radius   :: Radius
    } deriving Show

type Molecule = Map.Map ID Atom

instance Show Point where
    show (Point x y z) = unwords $ map show [x,y,z]

newMolecule :: Molecule
newMolecule = Map.fromList []

addAtom :: (ID, Atom) -> Molecule -> Molecule
addAtom (id, atom) molecule = Map.insert id atom molecule

-- | Установка рабочей области
workSpace :: [Point] -> [Point]
workSpace points = 
    let compX a@(Point x1 _ _) b@(Point x2 _ _) = compare x1 x2
        compY a@(Point _ y1 _) b@(Point _ y2 _) = compare y1 y2
        compZ a@(Point _ _ z1) b@(Point _ _ z2) = compare z1 z2
    in take 6 (List.sortBy compZ . List.sortBy compY . List.sortBy compX $ points)

-- | Функция определяет, находится ли атом
-- в рабочей области
isWorkSpace :: Atom -> [Point] -> Bool
isWorkSpace atom points =
    let r1:r2:r3:r4:r5:r6:r7:r8:_ = points
        (x, y, z) = (getX $ coordin atom, getY $ coordin atom, getZ $ coordin atom)
        between x a b = a < x && x < b
    in and [x `between` (getX r1) (getX r2), y `between` (getY r1) (getY r3), z `between` (getZ r1) (getZ r5)]


readMoleculeFromFile :: FilePath -> Molecule
readMoleculeFromFile inp = 
    let s = lines $ unsafePerformIO $ readFile inp
        atoms = fmap readAtom s
    in foldr addAtom Map.empty atoms 

readAtom :: String -> (ID, Atom)
readAtom str = 
    let a0:a1:a2:a3:a4:a5:a6:a7:a8:_ = words str
        id      = ID      $ read a0
        name    = Name    a1
        resname = ResName a2
        resseq  = ResSeq  $ read a3
        coordin = Point (read a4) (read a5) (read a6)
        element = Element a7
        radius  = Radius $ read a8
        atom    = Atom name resname resseq coordin element radius
    in  (id, atom)

reWriteFile :: FilePath -> Molecule
reWriteFile inf =  undefined

-- | Основная функцию.
-- Принцип работы.
-- На вход функции:
-- 1. Атом, который будем крепить
-- 2. Атом, к которому будет кропить
-- 3. Длина связи между атомами
-- После задания параметров выбирается рабочая область с центром в атоме, к которому крепим.
-- 
upgradeMolecule :: Atom -> Atom -> Double -> Molecule -> Molecule
upgradeMolecule sAtom eAtom dist molecule =
    let (sx, sy, sz) = (getX $ coordin sAtom, getY $ coordin sAtom, getZ $ coordin sAtom)
        r1' = Point (sx - 2 * dist) (sy - 2 * dist) (sz - 2 * dist) 
        r2' = Point (sx + 2 * dist) (sy - 2 * dist) (sz - 2 * dist) 
        r3' = Point (sx - 2 * dist) (sy + 2 * dist) (sz - 2 * dist) 
        r4' = Point (sx + 2 * dist) (sy + 2 * dist) (sz - 2 * dist) 
        r5' = Point (sx - 2 * dist) (sy - 2 * dist) (sz + 2 * dist) 
        r6' = Point (sx + 2 * dist) (sy - 2 * dist) (sz + 2 * dist)
        r7' = Point (sx - 2 * dist) (sy + 2 * dist) (sz + 2 * dist)
        r8' = Point (sx + 2 * dist) (sy + 2 * dist) (sz + 2 * dist)
        [r1,r2,r3,r4,r5,r6,r7,r8] = workSpace [r1',r2',r3',r4',r5',r6',r7',r8']
        molecule' = 
    in Map.empty

--reWriteFile :: FilePath -> Molecule
-- reWriteFile inf = 
--     let s = unsafePerformIO $ readFile inf 
--         newtxt = [ unwords $ a1:a2:a3:a4:a5:a6:a7:a8:a9:a10:a11:[] | 
--                     line <- lines s, 
--                     let a1:a2:a3:a4:a5:a6:a7:a8:a9:_:_:_:a10:[] = words line,
--                     let a11 = case a10 of "C" -> "0.356359487256"
--                                           "H" -> "0.040001352445"
--                                           "N" -> "0.329632525712"
--                                           "O" -> "0.302905564168"
--                                           "S" -> "0.356359487256"]
--     in writeFile "outinf" (unlines newtxt)


-- -- findRings :: Molecule -> [[Index]]
-- -- findRings molecule = 
-- --     let index = Map.keys  $ getBonds molecule
-- --         bonds = Map.elems $ getBonds molecule
-- --         filterEnum xs x
-- --             | (head xs == x) && (length xs > 2) = [xs ++ [x]]
-- --             | length xs > 7 = []
-- --             | x `elem` xs = []
-- --             | otherwise = createEnum (xs ++ [x]) $ bonds !! (x - 1)
-- --         createEnum xs1 xs2 = concat $ (\x -> filterEnum xs1 x) <$> xs2
-- --     in concat $ (\x -> createEnum [x] $ bonds !! (x - 1)) <$> index

-- addAtom :: Atom -> Molecule -> Molecule
-- addAtom atom molecule =
--     let index = 1 + (if Map.null serial then 0 else fst $ Map.findMax chain) 
--     in  molecule { getChain = Map.insert index atom chain }
    
-- addBonds :: Index -> [Index] -> Molecule -> Molecule
-- addBonds index conn molecule = 
--     let bonds = getBonds molecule
--         chain = getChain molecule 
--         conn' = filter (\x -> elem x $ Map.keys chain) conn 
--     in molecule { getBonds = Map.insert index conn' bonds}

-- deleteAtom :: Index -> Molecule -> Molecule
-- deleteAtom index molecule = 
--     let chain = getChain molecule 
--         bonds = getBonds molecule
--     in molecule { getChain = Map.delete index chain
--                 , getBonds = Map.map (\xs -> List.delete index xs) . Map.delete index $ bonds
--                 }