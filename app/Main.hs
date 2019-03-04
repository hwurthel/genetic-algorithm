module Main where

import Evolution
import Protein
import IO

main :: IO ()
main = print 1

-- Схема шага эволюции:
-- Генерация популяции -> Кроссинговер -> Мутации -> Селекция
-- Селекции также предшествует операция ранжирования. 