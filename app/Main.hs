module Main where

import Evolution
import Protein
import IO

main :: IO ()
main = print 1

-- Схема шага эволюции:
-- Кроссинговер -> Мутации -> Селекция