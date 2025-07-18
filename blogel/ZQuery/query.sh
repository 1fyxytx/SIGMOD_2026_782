rm -f run
make
OMP_NUM_THREADS=1 mpiexec -n 10 -quiet -cf /home/zyy/blogel/ZQuery/conf /home/zyy/blogel/ZQuery/run 