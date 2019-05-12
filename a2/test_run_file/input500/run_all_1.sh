make clean
make

for i in 1 2 4 8 16
do
 sbatch ./qsort_$i.sh
done

