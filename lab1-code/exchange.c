/**********************************************************************
 * Point-to-point communication using MPI
 *
 **********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  int rank, size;
  double a, b;
  int ierr;
  MPI_Status status;
  MPI_Request s_request, r_request;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */
  
  a = 100.0 + (double) rank;  /* Different a on different processors */

  if(rank == 0){
    MPI_Isend(&a, 1, MPI_DOUBLE, 1, 111, MPI_COMM_WORLD, &s_request);
    MPI_Irecv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &r_request);
    MPI_Wait(&s_request, &status);
    MPI_Wait(&r_request, &status);

    printf("Processor 0 got %f from processor 2\n", b);
    

  }
  else if (rank ==1){    
    MPI_Irecv(&b, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &r_request);
    MPI_Isend(&a, 1, MPI_DOUBLE, 0, 222, MPI_COMM_WORLD, &s_request);
    MPI_Wait(&s_request, &status);
    MPI_Wait(&r_request, &status);
        printf("Processor 1 got %f from processor 2\n", b);
    

  }



  
  /* Exchange variable a, notice the send-recv order */
/*  if (rank == 0) {
    MPI_Send(&a, 1, MPI_DOUBLE, 1, 111, MPI_COMM_WORLD);
    MPI_Recv(&b, 1, MPI_DOUBLE, 2, 333, MPI_COMM_WORLD, &status);
    printf("Processor 0 got %f from processor 2\n", b);
  } else if (rank==1) {
    MPI_Recv(&b, 1, MPI_DOUBLE, 0, 111, MPI_COMM_WORLD, &status);
    MPI_Send(&a, 1, MPI_DOUBLE, 2, 222, MPI_COMM_WORLD);
    printf("Processor 1 got %f from processor 0\n", b);
  }
  else if (rank == 2){
  	MPI_Recv(&b, 1, MPI_DOUBLE, 1, 222, MPI_COMM_WORLD, &status);
  	MPI_Send(&a, 1, MPI_DOUBLE, 0, 333, MPI_COMM_WORLD);
  	printf("Processor 2 got %f from processor 1\n", b);
  	
  }
  */
 
 //for(int i=0; i<size;i++){

  	// if(rank==0){
  	//     MPI_Send(&a, 1, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
   //  	MPI_Recv(&b, 1, MPI_DOUBLE, size-1, size-1, MPI_COMM_WORLD, &status);
   //  	printf("Processor %d got %f from processor %d\n", rank, b, size-1);
  	// }
  	// else if (rank == size-1){
  
	  // 	MPI_Recv(&b, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &status);
	  // 	MPI_Send(&a, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);
	  // 	printf("Processor %d got %f from processor %d\n", rank, b, rank-1); 		
  	// }
  	// else{
  	//   	MPI_Recv(&b, 1, MPI_DOUBLE, rank-1, rank-1, MPI_COMM_WORLD, &status);
  	// 	MPI_Send(&a, 1, MPI_DOUBLE, rank+1, rank, MPI_COMM_WORLD);
  	// 	printf("Processor %d got %f from processor %d\n", rank, b, rank-1);
  	
  	// }
  //	}



  printf("rank %d", size);


  MPI_Finalize(); 

  return 0;
}
