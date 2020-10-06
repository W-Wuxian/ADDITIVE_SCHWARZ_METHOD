 cbrange[0.8:2]
 dgrid3d,          32 ,          32
 hidden3d
 pm3d
 do for [i=           0 :           1 ]{splot 'SOL_NUMERIQUE/U.dat' index i u 1:2:3 with pm3d at b}
