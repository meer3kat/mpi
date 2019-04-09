/**********************************************************************
 * A simple "hello world" program for MPI/C
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
	int rank, size;

  	MPI_Init(&argc, &argv);               /* Initialize MPI               */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank); //get current process id
	MPI_Comm_size (MPI_COMM_WORLD, &size);
  //printf("Hello World!\n");             /* Print a message              */
	printf("hello world from process %d of total %d \n", rank, size);
  	MPI_Finalize();                       /* Shut down and clean up MPI   */

  	return 0;
}
