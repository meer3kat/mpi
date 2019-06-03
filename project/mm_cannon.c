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

	double timer = 0;
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
 	double* C_array = NULL;
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
        	C_array = (double*)malloc(n*n*sizeof(double *));
        }
        fclose(f);


        // mat_print(A_array, n);
        // mat_print(B_array, n);
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
 	
 // 	if(grid_rank == 0){
	//  	printf("rank: %d", grid_rank);
	//  	mat_print(myA, local_size);
	//  	mat_print(myB, local_size);
	//  	mat_print(myC, local_size);
	// }
	// printf("here i am after choping\n");

	// now we can strrt cannon's algorithms

	// mat_mult(myA, myB, myC, local_size);
	// mat_print(myC, local_size);

	int uprank, downrank, leftrank, rightrank;
	int shiftsource=0, shiftdest =0;
	// printf("shiftsource: %d\n", shiftsource);
	// printf("just initialize some neightbours\n");



	//initial alignment
	MPI_Cart_shift(comm_grid, 1, coordinates[0], &shiftsource, &shiftdest); //shift A
	// printf("shift A\n");
	// printf("shiftsource: %d", shiftsource);
	// printf("here i am after shiftg\n");
	MPI_Sendrecv_replace(myA, local_size*local_size, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_grid, &status);

	MPI_Cart_shift(comm_grid, 0, coordinates[1], &shiftsource, &shiftdest); //shift B
	MPI_Sendrecv_replace(myB, local_size*local_size, MPI_DOUBLE, shiftdest, 1, shiftsource, 1, comm_grid, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	// printf("finish initial shift A B	\n");

	MPI_Cart_shift(comm_grid, 1, 1, &rightrank, &leftrank);
	MPI_Cart_shift(comm_grid, 0, 1, &downrank, &uprank);

	for(int i=0; i<sqrt_size; i++){

		mat_mult(myA, myB, myC, local_size);
		// mat_print(myC, local_size);

		MPI_Sendrecv_replace(myA, local_size*local_size, MPI_DOUBLE, leftrank, 1, rightrank, 1, comm_grid, &status);

		MPI_Sendrecv_replace(myB, local_size*local_size, MPI_DOUBLE, uprank, 1, downrank, 1, comm_grid, &status);

	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(&(myC[0]), local_size*local_size, MPI_DOUBLE, C_array, sendcounts, displs, blocktype, 0, comm_grid);

	MPI_Barrier(MPI_COMM_WORLD);
	if(grid_rank==0){
		timer = MPI_Wtime() - timer;
		printf("%lf\n", timer);
		FILE * fp;
		fp = fopen("cannon_output.txt","a");
		fprintf (fp, "%d, %.8f, %d \n", n, timer, size);
		fclose(fp);
		mat_print(C_array, n);
	}



	MPI_Finalize ();
	return 0;
}

