set term png  

set pm3d map
set palette color



set output 'output,N=256.png'
N=256
set dgrid3d N+1,N+1

TOP=0.90
DY = 0.23

set key font ",10"

set autoscale
  
  filename="SOR_N=".N.".txt"
  set title "Model Problem Solution, N=".N
  set xlabel "x"
  set ylabel "y"
  splot filename using 1:2:3 
  unset title
  unset xlabel
  unset ylabel
  
unset output
