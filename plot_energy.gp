# plot_energy.gp

set datafile separator ','
set terminal pngcairo size 1200,800
set output "results/energy_leapfrog.png"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set grid
plot "results/Leapfrog_energy.csv" using 1:4 with lines title "Totale"

set output "results/energy_euler.png"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set grid
plot "results/Euler_energy.csv" using 1:4 with lines title "Totale"

set output "results/energy_BH.png"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set grid
plot "results/Barnes_Hut_energy.csv" using 1:4 with lines title "Totale"
