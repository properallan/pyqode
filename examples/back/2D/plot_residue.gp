set key autotitle columnhead 
set datafile separator ','
plot "history.csv" usin 1:2 with lines, '' usin 1:3 with lines, '' usin 1:4 with lines, '' usin 1:5 with lines, '' usin 1:6 with lines 
#pause
#reread