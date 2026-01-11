set title '2D Conduction in a Box'
set grid
set xlabel 'y-coord'
set ylabel 'Temperature'
p 'temp_yline.dat' u 1:($2-273) w l ti 'Present: x = 0.0833','x_pt0833.ref' u 1:2 ls 1 ti 'Reference: x = 0.0833'
rep 'temp_yline.dat' u 1:($3-273) w l ti 'Present: x = 0.5','x_pt5.ref' u 1:2 ls 3 ti 'Reference: x = 0.5'
