  set   autoscale                        # scale axes automatically
  unset log                              # remove any log-scaling
  unset label                            # remove any previous labels
  set xtic auto                          # set xtics automatically
  set ytic auto                          # set ytics automatically

  set title "CO"

  set xlabel "axial coordinate [cm]"
  set ylabel "CO mass fraction"
  
  set xr [0.:2]

set terminal png size 800,600 
set output 'CO.png'

  plot    "gri30/Output/Solution.final.out" using 2:89 title 'GRI3.0' with lines , \
		  "vc-co/Output/Solution.final.out" using 2:38  title 'VC' with lines 
 