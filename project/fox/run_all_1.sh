make clean
make

for i in 1 4 9 16 36 64 81 144 256
do
 sbatch ./fox_$i.sh
done

