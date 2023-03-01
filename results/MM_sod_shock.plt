reset
#set format y "10^{%L}
#set format x "10^{%L}
set terminal postscript  enhanced color "Helvetica,22"
set xlabel "x position (cm)" font "Helvetica,22"
set xtics font "Helvetica,20"
set ytics font "Helvetica,20"
#set datafile separator ","
#
set xrange[0:1]
set xtics 0,.1,1
#set yrange [0:1.05]
#
#set logscale y
#set logscale x

set grid
set style line 1 lt 1 lw 4 ps 1 pi 30
set style line 2 lt 3 lw 4 ps 1 pi 30
set style line 3 lt 4 lw 4 ps 1 pi 30
set style line 4 lt 5 lw 4 ps 1 pi 30
set style line 5 lt 6 lw 4 ps 1 pi 30
set style line 6 lt 7 lw 4 ps 1 pi 30

set tmargin 0.1
set rmargin 0.1
set bmargin 0.5
set lmargin 5

set key bottom center
set title "Moving Mesh Velocity"
set ylabel "Velocity " font "Helvetica,22"
set output "mm_shock_vel.eps"
plot "reimann.dat" u 1:4 w l ls 1 title "Analytical Velocity", "MM_sod_so.dat" u 1:3 w l ls 3 title "Second Order Velocity"

set key right top
set title "Moving Mesh Pressure"
set ylabel "Pressure" font "Helvetica,22"
set output "mm_shock_pres.eps"
plot "reimann.dat" u 1:3 w l ls 1 title "Analytical Pressure", "MM_sod_so.dat" u 1:5 w l ls 3 title "Second Order Pressure"

set key right top
set title "Moving Mesh Density"
set ylabel "Density" font "Helvetica,22"
set output "mm_shock_dens.eps"
plot "reimann.dat" u 1:2 w l ls 1 title "Analytical Density", "MM_sod_so.dat" u 1:2 w l ls 3 title "Second Order Density"

set title "Moving Mesh Volume"
set ylabel "Volume" font "Helvetica,22"
set output "mm_shock_vol.eps"
plot "MM_sod_so.dat" u 1:6 w l ls 3
