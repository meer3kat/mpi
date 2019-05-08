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
	int* arr = (int*)malloc(sizeof(int)*size);
	int i= rank;
	MPI_Allgather(&i, 1, MPI_INT, arr, 1, MPI_INT,MPI_COMM_WORLD);
	// if(rank==0){
	// 	for(int j=0; j<size;j++){
	// 		printf("%d ", arr[j]);
	// 	}
	// 	printf("\n");
	// }
	// if(rank==size-1){
	// 	for(int j=0; j<size;j++){
	// 		printf("%d ", arr[j]);
	// 	}
	// 	printf("\n");
	// }
	for(int j=0; j<size;j++){
			printf("%d ", arr[j]);
		}
		printf("\n");


	// Program arguments

	MPI_Finalize(); /* Shut down and clean up MPI */	
	return 0;

}
