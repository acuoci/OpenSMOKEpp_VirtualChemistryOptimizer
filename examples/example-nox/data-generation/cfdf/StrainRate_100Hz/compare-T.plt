  set   autoscale                        # scale axes automatically
  unset log                              # remove any log-scaling
  unset label                            # remove any previous labels
  set xtic auto                          # set xtics automatically
  set ytic auto                          # set ytics automatically

  set title "Temperature"

  set xlabel "axial coordinate [cm]"
  set ylabel "Temperature [K]"
  
  set xr [0.:2.0]
# set yr [0:325]

set terminal png size 800,600 
set output 'T.png'

  plot    "gri30/Output/Solution.final.out" using 2:3 title 'GRI3.0' with lines , \
		  "vc/Output/Solution.final.out" using 2:3 title 'VC' with lines 
 