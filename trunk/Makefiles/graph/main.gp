plot "gnuplot.dat" using 1:2 with lines title 'Best Evaluation', \
     "gnuplot.dat" using 1:3 with lines title 'Mean Evaluation', \
     "gnuplot.dat" using 1:4 with lines title 'Worst Evaluation'
load "loop.gp"
