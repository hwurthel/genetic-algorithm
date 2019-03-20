module Utils 
    (
        fst4, snd4,
        trd4, fth4
    )
where 

fst4 (x,_,_,_) = x
snd4 (_,x,_,_) = x
trd4 (_,_,x,_) = x
fth4 (_,_,_,x) = x