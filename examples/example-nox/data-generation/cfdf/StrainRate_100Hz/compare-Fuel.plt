# Gnuplot script file for plotting data in file "force.dat"

  set   autoscale                        # scale axes automatically
  unset log                              # remove any log-scaling
  unset label                            # remove any previous labels
  set xtic auto                          # set xtics automatically
  set ytic auto                          # set ytics automatically

  set title "Fuel"

  set xlabel "axial coordinate [cm]"
  set ylabel "FUEL mass fraction"
  
  set xr [0.:2]
#  set yr [0:325]

set terminal png size 800,600 
set output 'CH4.png'

  plot    "gri30/Output/Solution.final.out" using 2:88 title 'GRI3.0' with lines , \
		  "vc/Output/Solution.final.out" using 2:27 title 'VC' with lines 
 