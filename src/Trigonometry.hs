{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module Trigonometry where

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