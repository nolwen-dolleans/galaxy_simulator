# plot_energy.gp

set datafile separator ','
set terminal pngcairo size 1200,800
set output "results/energy_leapfrog.png"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set grid
plot "results/Leapfrog_energy.csv" using 1:2 with lines title "Cinétique", \
     "results/Leapfrog_energy.csv" using 1:3 with lines title "Potentielle", \
     "results/Leapfrog_energy.csv" using 1:4 with lines title "Totale"
