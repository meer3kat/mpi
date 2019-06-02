#include <math.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

// Print matrix
int read_file(char *name, double*** A, double*** B){
//return the number of number read in the file and pointer *pp will point to the first element
    FILE* f;
    f = fopen(name, "r");
    if(f){
        fseek(f, 0, SEEK_END);
        fseek(f, 0, SEEK_SET);
        double** mA = NULL;
        double** mB = NULL;
        int n;
        fscanf(f,"%d ",&n);
        mA = (double**)malloc(n*sizeof(double*));
        mB = (double**)malloc(n*sizeof(double*));

        for(int i=0;i<n;i++){
        	mA[i] = (double*)malloc(n*sizeof(double));
        	for(int j=0; j<n; j++){
        		fscanf(f, "%lf ", &mA[i][j]);
        	}
        }
        *A = mA;


        for(int i=0;i<n;i++){
        	mB[i] = (double*)malloc(n*sizeof(double));
        	for(int j=0; j<n; j++){
        		fscanf(f, "%lf ", &mB[i][j]);
        	}
        }
        *B = mB;       

        fclose(f);
        return n;
    }
    else{
    	printf("error with open your input file.\n");
    	return 0;
    }
}



void print_matrix(double **A, int n){
	printf ("\n");
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			printf("% 10.2f ", A[i][j]);
		}
		printf("; \n");
	}
	printf ("\n");

}

void mat_print (double *A, int sz)
{
  printf ("\n");
  for (int i = 0; i < sz; i++)
    {
      //printf ("|");
      for (int j = 0; j < sz; j++)
        printf ("% 10.2f", A[i * sz + j]);
      printf (";\n");
    }
  printf ("\n");
}

void mat_mult(double *A, double *B, double *C, int sz)
{
  for (int i = 0; i < sz; i++)
    {
      for (int j = 0; j < sz; j++)
        {
          for (int k = 0; k < sz; k++)
            {
              C[i * sz + j] += A[i * sz + k] * B[k * sz + j];
            }
        }
    }
}


int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv); //initialize 
	int rank, size, temp_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get my number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processors

	MPI_Status status;
	MPI_Request req_send, req_recv;
	double timer;
	// Program arguments
	if(argc != 3){
		printf("please enter 1 input; to run ./qsort inputfile outputfile methods\n.");
		return -1;
	}
	// printf("input file: %s\n",argv[0]);		
	// printf("pivot strategy %s\n", argv[3]);
	char* input_file = argv[1];
	char* output_file = argv[2];
	// printf("%s\n",input_file);
	double **A = NULL;
	double **B = NULL;
	double **C = NULL;
	int n = 0;
	// if(rank == 0){
	// 	n = read_file(input_file,&A, &B);

	// 	// print_matrix(A,n);
	// 	// print_matrix(B,n);
	// }
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD); //let all rank get n


	int sqrt_size = sqrt(size);

	if(sqrt_size*sqrt_size != size){
		if(rank == 0) {
			printf("need to run mpi on number of processor that has int sqrt\n");
		}
		MPI_Abort(MPI_COMM_WORLD, -1);
	}

	int dimensions[2], periods[2], coordinates[2], remain_dims[2];
	dimensions[0] = sqrt_size;
	dimensions[1] = sqrt_size;
	periods[0] = 1;
	periods[1] = 1;
	int ndim = 2;


	// create 2d cartesian grid
  	MPI_Comm comm_grid;
  	int grid_rank;
  	MPI_Dims_create(size, ndim, dimensions);
	MPI_Cart_create(MPI_COMM_WORLD, ndim, dimensions, periods, 1, &comm_grid);
  	MPI_Comm_rank(comm_grid, &grid_rank);
 	MPI_Cart_coords(comm_grid, grid_rank, ndim, coordinates);

 	// create row communicattor
 	MPI_Comm comm_row;
 	int row_rank;
 	MPI_Comm_split(comm_grid, coordinates[0], coordinates[1], &comm_row);
 	MPI_Comm_rank(comm_row, &row_rank);

 	// create column commniator

 	MPI_Comm comm_col;
 	int col_rank;
	MPI_Comm_split(comm_grid, coordinates[1], coordinates[0], &comm_col);
 	MPI_Comm_rank(comm_col, &col_rank);

 	double* A_array = NULL;
 	double* B_array = NULL;
 	if(grid_rank == 0){
 		FILE* f;
    	f = fopen(input_file, "r");
    	if(f){
	        fseek(f, 0, SEEK_END);
        	fseek(f, 0, SEEK_SET);
        	fscanf(f, "%d ", &n);
        	A_array = (double*)malloc(n*n*sizeof(double *));
        	for(int i=0; i<n*n; i++){
        		fscanf(f, "%lf ", &A_array[i]);
        	}
        	B_array = (double*)malloc(n*n*sizeof(double *));
        	for(int i=0; i<n*n; i++){
        		fscanf(f, "%lf ", &B_array[i]);
        	}
        }
        fclose(f);
        mat_print(A_array, n);
        mat_print(B_array, n);
 	}

 	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

 	int local_size;
 	local_size = n/sqrt_size;

 	//block size on each process
 	double *myA = NULL;
 	double *myB = NULL;
 	double *myC = NULL;

 	myA = (double*)malloc(local_size * local_size * sizeof(double));
 	myB = (double*)malloc(local_size * local_size * sizeof(double));
 	myC = (double*)malloc(local_size * local_size * sizeof(double));
 	memset(myC, 0, local_size * local_size * sizeof(double));

 	double* mytempA = NULL;
 	double* mytempB = NULL;

 	mytempA = (double*)malloc(local_size * local_size * sizeof(double));
 	mytempB = (double*)malloc(local_size * local_size * sizeof(double));

 	//create block type
 	int global_size[2] = {n, n};  //nxn matrix A and B
 	int block_size[2] = {local_size, local_size};  //local matrix block
 	int starts[2] = {0, 0}; //starting position of first block
 	MPI_Datatype type, blocktype;
 	MPI_Type_create_subarray(2, global_size, block_size, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
 	MPI_Type_create_resized(type, 0, local_size * sizeof(double), &blocktype);
 	MPI_Type_commit(&blocktype);

 	int sendcounts[sqrt_size * sqrt_size];
 	int displs[sqrt_size * sqrt_size];
 	int pos[2];
 	int block_start;
 	if(grid_rank == 0){
 		timer = MPI_Wtime();
 		for(int i=0; i<sqrt_size*sqrt_size; i++){
 			sendcounts[i] = 1;
 		}

 		for(int i=0; i<sqrt_size; i++){
 			for(int j=0; j<sqrt_size; j++){
 				pos[0] = i;
 				pos[1] = j;
 				MPI_Cart_rank(comm_grid, pos, &temp_rank);
 				block_start = i * local_size * n + j * local_size;
 				displs[temp_rank] = block_start/local_size;
 			}
 		}
 	}

 	MPI_Scatterv(&(A_array[0]), sendcounts, displs, blocktype, &(myA[0]), local_size*local_size, MPI_DOUBLE, 0, comm_grid);
 	MPI_Scatterv(&(B_array[0]), sendcounts, displs, blocktype, &(myB[0]), local_size*local_size, MPI_DOUBLE, 0, comm_grid);
 	
 	if(grid_rank == 0){
	 	printf("rank: %d", grid_rank);
	 	mat_print(myA, local_size);
	 	mat_print(myB, local_size);
	 	mat_print(myC, local_size);
	}
	printf("here i am after choping\n");

	// now we can strrt cannon's algorithms

	// mat_mult(myA, myB, myC, local_size);
	// mat_print(myC, local_size);

	int uprank, downrank, leftrank, rightrank;
	int shiftsource=0, shiftdest =0;
	printf("shiftsource: %d\n", shiftsource);
	printf("just initialize some neightbours\n");



	//initial alignment
	MPI_Cart_shift(comm_grid, 1, coordinates[0], &shiftsource, &shiftdest); //shift A
	printf("shift A\n");
	printf("shiftsource: %d", shiftsource);
	printf("here i am after shiftg\n");
	MPI_Sendrecv_replace(myA, local_size*local_size, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_grid, &status);

	MPI_Cart_shift(comm_grid, 0, coordinates[1], &shiftsource, &shiftdest); //shift B
	MPI_Sendrecv_replace(myB, local_size*local_size, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_grid, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	printf("finish initial shift A B	\n");

	MPI_Cart_shift(comm_grid, 1, 1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_grid, 0, 1, &downrank, &uprank);

	for(int i=0; i<sqrt_size; i++){

		mat_mult(myA, myB, myC, local_size);
		// mat_print(myC, local_size);

		MPI_Sendrecv_replace(myA, local_size*local_size, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_grid, &status);

		MPI_Sendrecv_replace(myB, local_size*local_size, MPI_DOUBLE, uprank, 1, downrank, 1, comm_grid, &status);

	}




	// for(cannon_block_cycle=0; cannon_block_cycle < sqrt_size; cannon_block_cycle++){
	// 	for(C_index=0, A_row = 0; A_row < local_size; A_row++){
	// 		for(B_column=0; B_column<local_size; B_column++, C_index++){
	// 			for(A_column=0; A_column < local_size; A_column++){
	// 				myC[C_index] += myA[A_row * local_size + A_column] * myB[A_column*local_size + B_column]; 
	// 			}
	// 		}
	// 	}

	// 	//rotate blocks horizontally
	// 	MPI_Sendrecv_replace(myA, local_size*local_size, MPI_DOUBLE, (coordinates[1]+sqrt_size-1)%sqrt_size, 0, (coordinates[1]+1)%sqrt_size, 0, comm_row, &status);
	// 	MPI_Sendrecv_replace(myB, local_size*local_size, MPI_DOUBLE, (coordinates[0]+sqrt_size-1)%sqrt_size, 0, (coordinates[0]+1)%sqrt_size, 0, comm_col, &status);
	// }

	MPI_Barrier(MPI_COMM_WORLD);
 	double* C_array = NULL;
 	C_array = (double*)malloc(n*n*sizeof(double));
	MPI_Gatherv(&(myC[0]), local_size*local_size, MPI_DOUBLE, C_array, sendcounts, displs, blocktype, 0, comm_grid);

	MPI_Barrier(MPI_COMM_WORLD);

	mat_print(C_array, n);






 	// // local arrays
 	// int A_local_block_rows,A_local_block_columns,A_local_size;
 	// double *A_local_block = NULL;
 	// A_local_block_rows = n/sqrt_size;
 	// A_local_block_columns = n/sqrt_size;
 	// A_local_size = A_local_block_rows * A_local_block_columns;
 	// A_local_block = (double*)malloc(A_local_size * sizeof(double));


 	// int B_local_block_rows,B_local_block_columns,B_local_size;
 	// double *B_local_block = NULL;
 	// B_local_block_rows = n/sqrt_size;
 	// B_local_block_columns = n/sqrt_size;
 	// B_local_size = B_local_block_rows * B_local_block_columns;
 	// B_local_block = (double*)malloc(B_local_size * sizeof(double));

 	// double *C_local_block = NULL;
 	// C_local_block = (double*)malloc(A_local_block_rows * B_local_block_columns * sizeof(double));
 	// for(int i=0; i<A_local_block_rows*B_local_block_columns; i++){
 	// 	C_local_block[i]=0;
 	// }

 	// double *A_array = NULL;
 	// double *B_array = NULL;
 	// double *C_array = NULL;
 	// int global_size[2] = {n, n};
 	// int local_size[2] = {n/sqrt_size, n/sqrt_size};
 	// int starts[2] = {0, 0};
 	// MPI_Datatype type, blocktype;
  // 	MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
  // 	MPI_Type_create_resized (type, 0, size * sizeof (double), &blocktype);
  // 	MPI_Type_commit (&blocktype);
  // 	int sendcounts[sqrt_size * sqrt_size]; //number of blocks received by each processor;
  // 	int displs[sqrt_size * sqrt_size]; //displacements of the starting positions of blocks
  // 	int pos[2];
  // 	int block_start;

 	// if(rank == 0){

 	// 	FILE* f;
  //   	f = fopen(input_file, "r");
  //   	if(f){
	 //        fseek(f, 0, SEEK_END);
  //       	fseek(f, 0, SEEK_SET);
  //       	fscanf(f, "%d ", &n);
  //       	A_array = (double*)malloc(n*n*sizeof(double *));
  //       	for(int i=0; i<n*n; i++){
  //       		fscanf(f, "%lf ", &A_array[i]);
  //       	}
  //       	B_array = (double*)malloc(n*n*sizeof(double *));
  //       	for(int i=0; i<n*n; i++){
  //       		fscanf(f, "%lf ", &B_array[i]);
  //       	}
  //       }
  //       fclose(f);

  //       timer = MPI_Wtime();

  //       for(int i=0; i<sqrt_size*sqrt_size; i++){
  //       	sendcounts[i] = 1;
  //       }

  //       for(int i=0; i<sqrt_size; i++){
  //       	for(int j=0; j<sqrt_size; j++){
  //       		pos[0] = i;
  //       		pos[1] = j;
  //       		MPI_Cart_rank(comm_grid, pos, &temp_rank);
  //       		block_start = i*local_size * n + j * local_size;
  //       		displs[temp_rank] = block_start / local_size;
  //       	}
  //       }
    
 	// }

 	// double* myA = (double*)malloc(local_size*local_size*sizeof(double));
 	// double* myB = (double*)malloc(local_size*local_size*sizeof(double));
 	// double* myC = (double*)malloc(local_size*local_size*sizeof(double));
 	// memset(myC,0, local_size*local_size*sizeof(double));
 	
 	



 	// MPI_Scatterv(&(A[0]), sendcounts, displs, blocktype, &(myA[0]),local_size*local_size, MPI_DOUBLE, 0, comm_grid);
 	// MPI_Scatterv(&(B[0]), sendcounts, displs, blocktype, &(myB[0]),local_size*local_size, MPI_DOUBLE, 0, comm_grid);
 	

 	// if(rank == 0){
 	// 	A_array = (double*)malloc(n*n*sizeof(double));
 	// 	B_array = (double*)malloc(n*n*sizeof(double));
 	// 	C_array = (double*)malloc(n*n*sizeof(double));

 	// 	for(int i=0; i < sqrt_size; i++){
 	// 		for(int j=0; j<sqrt_size; j++){
 	// 			for(int row=0; row<A_local_block_rows; row++){
 	// 				for(int column=0; column<A_local_block_columns;column++){
 	// 					A_array[((i*sqrt_size+j) * A_local_size) + (row * A_local_block_columns) + column] = A[i*A_local_block_rows + row][j*A_local_block_columns + column];
 	// 					B_array[((i*sqrt_size+j) * B_local_size) + (row * B_local_block_columns) + column] = B[i*B_local_block_rows + row][j*B_local_block_columns + column];
 	// 				}
 	// 			}
 	// 		}
 	// 	}

		// C = (double **)malloc(n * sizeof(double *));
		// for(int i=0; i<n ;i++){
		// 	C[i] = (double *)malloc(n * sizeof(double));
		// }

 	// }
 	// //send a block to each
 	// if(rank == 0){
 	// 	for(int i=0; i<size; i++){
 	// 		MPI_Send((A_array + (i*A_local_size)), A_local_size, MPI_Double, i, 0, comm_grid);
 	// 		MPI_Send
 	// 	}

 	// }

 	







	MPI_Finalize ();


	return 0;
}








// void *mem_alloc (int);
// void mem_reset (double **);
// void mat_fill (double *, int, FILE *);
// void mat_mult (double *, double *, double *, int);
// void mat_print (double *, int);
// void mat_fprint (double *, int, char *);

// int main (int argc, char *argv[])
// {	
// 	//read file 





//   FILE *inFilePtr;
//   int p, sqrtp, N, size, temp_rank, block_start;
//   int grid_rank, row_rank, col_rank;
//   int bcast_root, source, dest;
//   double *myA, *myB, *myC, *mytempA, *mytempB;
//   double *A, *B, *C;
//   double timer;

//   int coords[2], pos[2], reorder = 1, ndim = 2;
//   int dims[2] = { 0, 0 };
//   int periods[2] = { 0, 0 };

//   MPI_Status status;
//   MPI_Request req_send, req_recv;

//   // Initialize MPI
//   MPI_Init (&argc, &argv);
//   MPI_Comm_size (MPI_COMM_WORLD, &p);

//   sqrtp = (int) sqrt ((double) p); //square-root of number of processors
//   dims[0] = sqrtp;
//   dims[1] = sqrtp;

//   // Create a 2D Cartesian topology
//   MPI_Comm comm_grid;
//   MPI_Dims_create (p, ndim, dims);
//   MPI_Cart_create (MPI_COMM_WORLD, ndim, dims, periods, reorder, &comm_grid);
//   MPI_Comm_rank (comm_grid, &grid_rank);
//   MPI_Cart_coords (comm_grid, grid_rank, ndim, coords);

//   // Create a communicator for each row
//   MPI_Comm comm_row;
//   MPI_Comm_split (comm_grid, coords[0], coords[1], &comm_row);
//   MPI_Comm_rank (comm_row, &row_rank);

//   // Create a communicator for each column
//   MPI_Comm comm_col;
//   MPI_Comm_split (comm_grid, coords[1], coords[0], &comm_col);
//   MPI_Comm_rank (comm_col, &col_rank);

//   MPI_Barrier (MPI_COMM_WORLD);

//   if (grid_rank == 0)
//     {
//         char *inFileName = argv[1];
//         inFilePtr = fopen(inFileName, "r");
//         fscanf(inFilePtr, "%d", &N);
//     }
//   /* Broadcast value of N to all the other nodes */
//   MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

  
//   size = N / sqrtp; //order of each block

//   // Allocate memory for the block matrices
//   myA = mem_alloc (size * size * sizeof (double));
//   myB = mem_alloc (size * size * sizeof (double));
//   myC = mem_alloc (size * size * sizeof (double));
//   memset (myC, 0, size * size * sizeof (double));

//   mytempA = mem_alloc (size * size * sizeof (double));
//   mytempB = mem_alloc (size * size * sizeof (double));

//   // Create block type
//   int sizes[2] = { N, N };      //global size of matrix
//   int subsizes[2] = { size, size };     //size of each block
//   int starts[2] = { 0, 0 };     //position of first block
//   MPI_Datatype type, blocktype;
//   MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
//   MPI_Type_create_resized (type, 0, size * sizeof (double), &blocktype);
//   MPI_Type_commit (&blocktype);

//   int sendcounts[sqrtp * sqrtp]; //number of blocks received by each processor;
//   int displs[sqrtp * sqrtp]; //displacements of the starting positions of blocks

//   if (grid_rank == 0)
//     {
//       // Allocate memory for the "full" matrices
//       A = mem_alloc (N * N * sizeof (double));
//       B = mem_alloc (N * N * sizeof (double));
//       C = mem_alloc (N * N * sizeof (double));

//       // Supplies current time as seed for random number generation
//       srand (time (NULL));

//       // Generate random matrices A and B
//       mat_fill (A, N, inFilePtr);
//       mat_fill (B, N, inFilePtr);
//       fclose(inFilePtr);
      
//       // Print A and B for debugging
//       //printf ("\nA=");
//       //mat_print (A, N);
//       //printf ("\nB=");
//       //mat_print (B, N);

//       timer = MPI_Wtime ();
//       for (int i = 0; i < sqrtp * sqrtp; i++)
//         sendcounts[i] = 1;      //each processor receives one block
//       //processor (i,j) receives block (i,j)
//       for (int i = 0; i < sqrtp; i++)
//         {
//           for (int j = 0; j < sqrtp; j++)
//             {
//               pos[0] = i;
//               pos[1] = j;
//               MPI_Cart_rank (comm_grid, pos, &temp_rank);
//               block_start = i * size * N + j * size;
//               // compute the starting point of each processors block
//               // in the global matrix, in block extents 
//               displs[temp_rank] = block_start / size;
//             }
//         }
//     }
//   // Scatter matrices A and B to all the processors
//   MPI_Scatterv (&(A[0]), sendcounts, displs, blocktype, &(myA[0]), size * size, MPI_DOUBLE, 0, comm_grid);
//   MPI_Scatterv (&(B[0]), sendcounts, displs, blocktype, &(myB[0]), size * size, MPI_DOUBLE, 0, comm_grid);

// 	// Delete matrices A and B in the root node, which are not needed any more
//   if (grid_rank == 0)
//     {
//       mem_reset (&A);
//       mem_reset (&B);
//     }

//   source = (col_rank + 1) % sqrtp; //source node when shifting B 
//   dest = (col_rank + sqrtp - 1) % sqrtp; //destination node when shifting B

//   MPI_Barrier (MPI_COMM_WORLD);

//   for (int k = 0; k < sqrtp; k++)
//     {
//       bcast_root = (col_rank + k) % sqrtp; //compute broadcasting node for each row
//       if (bcast_root == row_rank) //the broadcasting processor needs the right matrix too
//         memcpy (&mytempA[0], &myA[0], size * size * sizeof (double));
//       // Broadcast A to all the processors in comm_row
//       MPI_Bcast (&mytempA[0], size * size, MPI_DOUBLE, bcast_root, comm_row);

//       // Shift B using non-blocking communication
//       MPI_Irecv (&mytempB[0], size * size, MPI_DOUBLE, source, 666, comm_col, &req_recv);
//       MPI_Isend (&myB[0], size * size, MPI_DOUBLE, dest, 666, comm_col, &req_send);

//       // Multiply matrices
//       mat_mult (mytempA, myB, myC, size);
			
//       // Wait for receiving the "new" B matrix to arrive
//       MPI_Wait (&req_recv, &status);
//       MPI_Wait (&req_send, &status);

//       memcpy (&myB[0], &mytempB[0], size * size * sizeof (double));
//     }

//   MPI_Barrier (MPI_COMM_WORLD);

// 	//Build the resultant matrix C
//   MPI_Gatherv (&(myC[0]), size * size, MPI_DOUBLE, C, sendcounts, displs, blocktype, 0, comm_grid);

//   MPI_Barrier (MPI_COMM_WORLD);
//   if (grid_rank == 0)
//     {
//       timer = MPI_Wtime () - timer;
//       printf ("%lf\n",timer);
//     }
//   mem_reset (&myA);
//   mem_reset (&myB);
//   mem_reset (&myC);
//   mem_reset (&mytempA);
//   mem_reset (&mytempB);
// 	///*
//   if (grid_rank == 0)
//     {
//       //printf ("\nC=");
//       //mat_print (C, N);
//       mat_fprint (C, N, argv[2]);
//       mem_reset (&C);
//     }
// 	//*/
//   MPI_Type_free (&blocktype);
//   MPI_Comm_free (&comm_row);
//   MPI_Comm_free (&comm_col);
//   MPI_Comm_free (&comm_grid);

//   MPI_Finalize ();
//   return 0;
// }

// // Dynamically allocate memory for an array
// void *
// mem_alloc (int n)
// {
//   void *ptr_mem = malloc (n);
//   if (ptr_mem == NULL)
//     {
//       printf ("\nERROR! Not enough memory available.");
//       MPI_Finalize ();
//       exit (EXIT_FAILURE);
//     }

//   return ptr_mem;
// }

// // Free memory of dynamically allocated array
// void
// mem_reset (double **ptr)
// {
//   if (ptr)
//     {
//       free (*ptr);
//       *ptr = NULL;
//     }
// }

// // Fill matrix with a random double
// void
// mat_fill (double *A, int sz, FILE* fptr)
// {
//   for (int i = 0; i < sz * sz; i++)
//     fscanf(fptr, "%lf", &A[i]);
// }

// // Multiply matrix blocks
// void
// mat_mult (double *A, double *B, double *C, int sz)
// {
//   for (int i = 0; i < sz; i++)
//     {
//       for (int j = 0; j < sz; j++)
//         {
//           for (int k = 0; k < sz; k++)
//             {
//               C[i * sz + j] += A[i * sz + k] * B[k * sz + j];
//             }
//         }
//     }
// }

// void
// // Print matrix
// mat_print (double *A, int sz)
// {
//   printf ("\n");
//   for (int i = 0; i < sz; i++)
//     {
      //printf ("|");
//       for (int j = 0; j < sz; j++)
//         printf ("% 10.2f", A[i * sz + j]);
//       printf (";\n");
//     }
//   printf ("\n");
// }

// void
// // Print matrix to file
// mat_fprint (double *A, int sz, char *outFileName)
// {
//   FILE *outFilePtr;
//   outFilePtr = fopen(outFileName, "w");

//   for (int i = 0; i < sz; i++)
//     {
//       for (int j = 0; j < sz; j++)
//         fprintf (outFilePtr, "%.6f ", A[i * sz + j]);
//       fprintf (outFilePtr,"\n");
//     }
//   fclose(outFilePtr);
// }
