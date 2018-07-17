  set   autoscale                        # scale axes automatically
  unset log                              # remove any log-scaling
  unset label                            # remove any previous labels
  set xtic auto                          # set xtics automatically
  set ytic auto                          # set ytics automatically

  set title "NO"

  set xlabel "axial coordinate [cm]"
  set ylabel "NO mass fraction"
  
  set xr [0.:2]

set terminal png size 800,600 
set output 'NO.png'

  plot    "gri30/Output/Solution.final.out" using 2:110 title 'GRI3.0' with lines , \
		  "vc-no-20180605/Output/Solution.final.out" using 2:39  title 'VC-20180605' with lines, \
		  "vc-no-20180619/Output/Solution.final.out" using 2:39  title 'VC-20180619' with lines
 