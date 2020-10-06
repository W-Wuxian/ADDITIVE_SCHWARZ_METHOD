 set cbrange[0:1]
 set dgrid3d,          22 ,          22
 set hidden3d
 set pm3d
 do for [i=           1 :         200 ]{splot 'SOL_NUMERIQUE/U.dat' index i u 1:2:3 with pm3d at sb}
