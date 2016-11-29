#!/usr/bin/gnuplot

dx="0.100"
index=5

#
#	PLOT THE LEARNING CURVE (AVERAGE COST)
#
reset
set term qt 0 persist
p dir.'avg_cost_dx'.dx.'.dat' u 1:2 w l, '' u 1:3 w l lw 2


#
#	PLOT THE DENSITY AND THE VALUE FIELD
#
reset
set view map; set pm3d map; unset contours; set size ratio -1

set term qt 1
splot dir.'rho_value_dx'.dx.'.dat' i index u 1:2:3

set term qt 2
splot dir.'rho_value_dx'.dx.'.dat' i index u 1:2:4 



#
#	PLOT THE DRIFT FIELD
#
reset
set term qt 3
set size ratio -1
p dir.'drift_dx'.dx.'.dat' i index u 1:2:3:4 w vec 


#
#	PLOT TRAJECTORY
#
reset
# png
set terminal pngcairo size 350,350 enhanced font 'Verdana,10'
# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue
unset key
set size ratio -1
set xrange [-10:10]
set yrange [-10:10]

n=0
do for [ii=1:99] {
    n=n+1
    set output sprintf(dir.'trajs/path_%03.0f.png',n)
    p dir.'typlen_0.40/N_1.0e4/g_0.000/trajs/traj_0001.dat' every ::1::ii u 2:3 w l ls 1, \
      '' every ::ii::ii u 2:3 w p ls 1
}




pause -1
