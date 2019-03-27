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

-- | Синонимы типов
type Elem  = String
type Index = Int
--type Param a = Map.Map Index a

data Atom = Atom
    { getElement :: Elem
    , getCoord   :: Point
    }

data Molecule = Molecule
    { getChain :: Map.Map Index Atom
    , getBonds :: Map.Map Index [Index]
    }

instance Show Point where
    show (Point x y z) = unwords $ map show [x,y,z]
                         
instance Show Atom where
    show atom =    (show $ getElement atom) <> "\t" 
                <> (show $ getCoord   atom)

instance Show Molecule where
    show molecule =  (show $ Map.toList $ getChain molecule) <> "\n" 
                  <> (show $ Map.toList $ getBonds molecule)

newMolecule :: Molecule
newMolecule = Molecule (Map.fromList []) (Map.fromList [])

findRings :: Molecule -> [[Index]]
findRings molecule = 
    let index = Map.keys  $ getBonds molecule
        bonds = Map.elems $ getBonds molecule
        filterEnum xs x
            | (head xs == x) && (length xs > 2) = [xs ++ [x]]
            | length xs > 7 = []
            | x `elem` xs = []
            | otherwise = createEnum (xs ++ [x]) $ bonds !! (x - 1)
        createEnum xs1 xs2 = concat $ (\x -> filterEnum xs1 x) <$> xs2
    in concat $ (\x -> createEnum [x] $ bonds !! (x - 1)) <$> index

addAtom :: Atom -> Molecule -> Molecule
addAtom atom molecule =
    let chain = getChain molecule
        index = 1 + (if Map.null chain then 0 else fst $ Map.findMax chain) 
    in  molecule { getChain = Map.insert index atom chain }
    
addBonds :: Index -> [Index] -> Molecule -> Molecule
addBonds index conn molecule = 
    let bonds = getBonds molecule
        chain = getChain molecule 
        conn' = filter (\x -> elem x $ Map.keys chain) conn 
    in molecule { getBonds = Map.insert index conn' bonds}

deleteAtom :: Index -> Molecule -> Molecule
deleteAtom index molecule = 
    let chain = getChain molecule 
        bonds = getBonds molecule
    in molecule { getChain = Map.delete index chain
                , getBonds = Map.map (\xs -> List.delete index xs) . Map.delete index $ bonds
                }

readMolecule :: FilePath -> Molecule
readMolecule = undefined
    -- 1. 