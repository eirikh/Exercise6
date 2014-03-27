set terminal post eps
set output "graph.eps"
set title "2D plot of a" 
splot 'matrix_output.dat' binary format="%lf" array=127x127 dx=0.0078125 dy=0.0078125 using 1 with lines
