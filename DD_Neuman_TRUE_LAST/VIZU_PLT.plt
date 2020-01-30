 set cbrange[   9.6407132793526569E-005 :   6.3222777676939471E-002 ]
 set dgrid3d,         102 ,         102
 set hidden3d
 set pm3d
 splot 'SOL_NUMERIQUE/U.dat'  u 1:2:3 with pm3d at sb
