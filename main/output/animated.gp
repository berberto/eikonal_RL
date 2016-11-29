#!/bin/gnuplot

reset

# png
set terminal pngcairo size 640,640 enhanced font 'Verdana,10'

# color definitions
set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 ps 2 # --- blue

unset key
set border 4095
set pm3d map; set size ratio -1; unset contours
set xrange [-3:3]
set yrange [-3:3]


dir = '/media/DATI/Fisica/DATI-EIKONAL/typlen_0.50/N_4.0e4/g_0.0000/trajs/'
#dir = 'typlen_0.50/N_4.0e4/g_0.0000/trajs/'

tr = '0002'
file=sprintf('traj_%s.dat', tr)

n=0
do for [ii=1:1000] {
    n=n+1
    set output sprintf('%s%03.0f.png', dir, n)
    p dir.file every ::1::ii u 2:3 w l ls 1 lt 1,  dir.file every ::ii::ii u 2:3 w p ls 1 pt 2
}

# Create movie with mencoder
ENCODER = system('which mencoder');
if (strlen(ENCODER)==0) print '=== mencoder not found ==='; exit
CMD = sprintf('mencoder mf://%s/*.png -mf fps=25:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o traj.avi',dir)
system(CMD)

system(sprintf("rm %s/*.png", dir))
