module Main where

import Evolution
import Protein
import IO

main :: IO ()
main = print 1

-- Схема шага эволюции:
-- Кроссинговер -> Мутации -> Селекция
-- !!! Необходимо не забыть реализовать
-- !!! запись всех получаемых белков в процессе работы