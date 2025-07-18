rm -f run
make
OMP_NUM_THREADS=4 mpiexec -n 10 -quiet -cf /home/zyy/blogel/Zdistance/conf /home/zyy/blogel/Zdistance/run 