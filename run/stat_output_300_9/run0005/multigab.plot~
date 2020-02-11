set terminal pngcairo size 1200, 1600
set output 'stacking.png'
set lmargin at screen 0.08
set rmargin at screen 0.95
TOP=0.95
DY = 0.22

set multiplot
unset key
set xrange[0:2]
set xtics font ",24"
set ytics font ",24"
set offset 0,0,graph 0.05, graph 0.05

set xlabel '{/Symbol D}{/Symbol h}, Spatial Rapidity' offset 0,-2 enhanced font ",32"
set ylabel 'SS' offset -3 font ",36"
set ytics ("-500" -500, "" 0, "500" 500,"" 1000,"1500" 1500,"2000" 2000,"2500" 2500)
set tmargin at screen TOP-3*DY
set bmargin at screen TOP-4*DY
plot for[i=5:801:4] "gabfunctions.dat" u 1:i lc rgb 'red' w l, for[i=5:81:4] "posteriorgabfunctions.dat" u 1:i lc rgb 'blue' w l, "gabfunctions.dat" u 1:25 lw 3 lc rgb 'black' w l

set xtics format ''
unset xlabel
set ylabel 'US' offset -3 font ",36"
set ytics ("-500" -500, "" -400,"-300" -300,"" -200,"-100" -100,"" 0,"100" 100)
set tmargin at screen TOP-2*DY
set bmargin at screen TOP-3*DY
plot for[i=4:801:4] "gabfunctions.dat" u 1:i lc rgb 'red' w l, for[i=4:81:4] "posteriorgabfunctions.dat" u 1:i lc rgb 'blue' w l, "gabfunctions.dat" u 1:24 lw 5 lc rgb 'black' w l

set ylabel 'UD' offset -4 font ",36"
set ytics ("-1500" -1500, "" -1000, "-500" -500,"" 0)
set tmargin at screen TOP-1*DY
set bmargin at screen TOP-2*DY
plot for[i=3:801:4] "gabfunctions.dat" u 1:i lc rgb 'red' w l, for[i=3:81:4] "posteriorgabfunctions.dat" u 1:i lc rgb 'blue' w l, "gabfunctions.dat" u 1:23 lw 5 lc rgb 'black' w l

set ylabel 'UU' offset -4 font ",36"
set ytics ("-1000" -1000,"" 0,"1000" 1000,"" 2000,"3000" 3000,"" 4000,"5000" 5000)
set tmargin at screen TOP-0*DY
set bmargin at screen TOP-1*DY
set label "Prior" textcolor rgb "red" font ",36" at 1.5,3500
set label "Posterior" textcolor rgb "blue" font ",36" at 1.5,2750
set label "Full Calculation" textcolor rgb "black" font ",36" at 1.5,2000
plot for[i=2:801:4] "gabfunctions.dat" u 1:i lc rgb 'red' w l, for[i=2:81:4] "posteriorgabfunctions.dat" u 1:i lc rgb 'blue' w l, "gabfunctions.dat" u 1:22 lw 5 lc rgb 'black' w l

unset multiplot; set output