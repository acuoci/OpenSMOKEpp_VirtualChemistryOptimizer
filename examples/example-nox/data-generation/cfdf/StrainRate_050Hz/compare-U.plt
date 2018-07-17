# Gnuplot script file for plotting data in file "force.dat"

  set   autoscale                        # scale axes automatically
  unset log                              # remove any log-scaling
  unset label                            # remove any previous labels
  set xtic auto                          # set xtics automatically
  set ytic auto                          # set ytics automatically

  set title "Axial velocity"

  set xlabel "axial coordinate [cm]"
  set ylabel "axial velocity [cm/s]"
  
  set xr [0.:2]

set terminal png size 800,600 
set output 'U.png'

  plot    "gri30/Output/Solution.final.out" using 2:4 title 'GRI3.0' with lines , \
		  "vc/Output/Solution.final.out" using 2:4 title 'VC' with lines 
 