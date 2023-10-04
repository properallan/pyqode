set key autotitle columnhead 
set datafile separator ','
plot "history.csv" usin 1:7 with lines, '' usin 1:8 with lines, '' usin 1:9 with lines 
#plot "history.csv" usin 1:9 with lines 
pause
reread