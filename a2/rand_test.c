/*****************************************************
* Run by typing: "mpirun -np p ./qsort"    		 *
* p: Number of processors (square number)            *
*****************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>



int main(int argc, char *argv[]){	
	// set up
	MPI_Init(&argc, &argv); //initialize 
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get my number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processors
	unsigned int key;

	// Program arguments

	if (rank == 0){
		srand(time(NULL));   // Initialization, should only be called once.
		key = rand();      // Returns a pseudo-random integer between 0 and RAND_MAX.
	}

	MPI_Bcast(&key, 1, MPI_INT,0, MPI_COMM_WORLD);
	unsigned int loc_seed;
	loc_seed = key+rank;

	int loc_rand;
	loc_rand = rand_r(&loc_seed);
	printf("thread %d, rand: %d\n",rank, loc_rand);

	MPI_Finalize(); /* Shut down and clean up MPI */	
	return 0;

}
