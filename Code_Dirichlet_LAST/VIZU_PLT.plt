 set cbrange[   1.8844599181223784E-004 :   6.2469846285770958E-002 ]
 set dgrid3d,         102 ,          52
 set hidden3d
 set pm3d
 splot 'SOL_NUMERIQUE/U.dat'  u 1:2:3 with pm3d at sb
