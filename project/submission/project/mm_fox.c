#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void *mem_alloc (int);
void mem_reset (double **);
void mat_fill (double *, int, FILE *);
void mat_mult (double *, double *, double *, int);
void mat_print (double *, int);
void mat_fprint (double *, int, char *);

int
main (int argc, char *argv[])
{
  FILE *inFilePtr;
  int p, sqrtp, N, size, temp_rank, block_start;
  int grid_rank, row_rank, col_rank;
  int bcast_root, source, dest;
  double *myA, *myB, *myC, *mytempA, *mytempB;
  double *A, *B, *C;

  double timer = 0.0;
  double start_time = 0.0;
  double read_file_time = 0.0;
  double all_time = 0.0;
  double prep_time = 0.0;

  int coords[2], pos[2], reorder = 1, ndim = 2;
  int dims[2] = { 0, 0 };
  int periods[2] = { 0, 0 };

  MPI_Status status;
  MPI_Request req_send, req_recv;

  // Initialize MPI
  MPI_Init (&argc, &argv);
  MPI_Comm_size (MPI_COMM_WORLD, &p);

  sqrtp = (int) sqrt ((double) p); //square-root of number of processors
  dims[0] = sqrtp;
  dims[1] = sqrtp;

  // Create a 2D Cartesian topology
  MPI_Comm comm_grid;
  MPI_Dims_create (p, ndim, dims);
  MPI_Cart_create (MPI_COMM_WORLD, ndim, dims, periods, reorder, &comm_grid);
  MPI_Comm_rank (comm_grid, &grid_rank);
  MPI_Cart_coords (comm_grid, grid_rank, ndim, coords);

  // Create a communicator for each row
  MPI_Comm comm_row;
  MPI_Comm_split (comm_grid, coords[0], coords[1], &comm_row);
  MPI_Comm_rank (comm_row, &row_rank);

  // Create a communicator for each column
  MPI_Comm comm_col;
  MPI_Comm_split (comm_grid, coords[1], coords[0], &comm_col);
  MPI_Comm_rank (comm_col, &col_rank);

  MPI_Barrier (MPI_COMM_WORLD);

  if (grid_rank == 0)
    {
        char *inFileName = argv[1];
        inFilePtr = fopen(inFileName, "r");
        fscanf(inFilePtr, "%d", &N);
        start_time = MPI_Wtime();
    }
  /* Broadcast value of N to all the other nodes */
  MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

  
  size = N / sqrtp; //order of each block

  // Allocate memory for the block matrices
  myA = mem_alloc (size * size * sizeof (double));
  myB = mem_alloc (size * size * sizeof (double));
  myC = mem_alloc (size * size * sizeof (double));
  memset (myC, 0, size * size * sizeof (double));

  mytempA = mem_alloc (size * size * sizeof (double));
  mytempB = mem_alloc (size * size * sizeof (double));

  // Create block type
  int sizes[2] = { N, N };      //global size of matrix
  int subsizes[2] = { size, size };     //size of each block
  int starts[2] = { 0, 0 };     //position of first block
  MPI_Datatype type, blocktype;
  MPI_Type_create_subarray (2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &type);
  MPI_Type_create_resized (type, 0, size * sizeof (double), &blocktype);
  MPI_Type_commit (&blocktype);

  int sendcounts[sqrtp * sqrtp]; //number of blocks received by each processor;
  int displs[sqrtp * sqrtp]; //displacements of the starting positions of blocks

  if (grid_rank == 0)
    {
      // Allocate memory for the "full" matrices
      A = mem_alloc (N * N * sizeof (double));
      B = mem_alloc (N * N * sizeof (double));
      C = mem_alloc (N * N * sizeof (double));

      // Supplies current time as seed for random number generation
      srand (time (NULL));

      // Generate random matrices A and B
      mat_fill (A, N, inFilePtr);
      mat_fill (B, N, inFilePtr);
      fclose(inFilePtr);
      read_file_time = MPI_Wtime() - start_time;
      
      // Print A and B for debugging
      //printf ("\nA=");
      //mat_print (A, N);
      //printf ("\nB=");
      //mat_print (B, N);

      timer = MPI_Wtime ();
      for (int i = 0; i < sqrtp * sqrtp; i++)
        sendcounts[i] = 1;      //each processor receives one block
      //processor (i,j) receives block (i,j)
      for (int i = 0; i < sqrtp; i++)
        {
          for (int j = 0; j < sqrtp; j++)
            {
              pos[0] = i;
              pos[1] = j;
              MPI_Cart_rank (comm_grid, pos, &temp_rank);
              block_start = i * size * N + j * size;
              // compute the starting point of each processors block
              // in the global matrix, in block extents 
              displs[temp_rank] = block_start / size;
            }
        }
      prep_time = MPI_Wtime() - start_time - read_file_time;
    }
  // Scatter matrices A and B to all the processors
  MPI_Scatterv (&(A[0]), sendcounts, displs, blocktype, &(myA[0]), size * size, MPI_DOUBLE, 0, comm_grid);
  MPI_Scatterv (&(B[0]), sendcounts, displs, blocktype, &(myB[0]), size * size, MPI_DOUBLE, 0, comm_grid);

	// Delete matrices A and B in the root node, which are not needed any more
  if (grid_rank == 0)
    {
      mem_reset (&A);
      mem_reset (&B);
    }

  source = (col_rank + 1) % sqrtp; //source node when shifting B 
  dest = (col_rank + sqrtp - 1) % sqrtp; //destination node when shifting B

  MPI_Barrier (MPI_COMM_WORLD);

  for (int k = 0; k < sqrtp; k++)
    {
      bcast_root = (col_rank + k) % sqrtp; //compute broadcasting node for each row
      if (bcast_root == row_rank) //the broadcasting processor needs the right matrix too
        memcpy (&mytempA[0], &myA[0], size * size * sizeof (double));
      // Broadcast A to all the processors in comm_row
      MPI_Bcast (&mytempA[0], size * size, MPI_DOUBLE, bcast_root, comm_row);

      // Shift B using non-blocking communication
      MPI_Irecv (&mytempB[0], size * size, MPI_DOUBLE, source, 666, comm_col, &req_recv);
      MPI_Isend (&myB[0], size * size, MPI_DOUBLE, dest, 666, comm_col, &req_send);

      // Multiply matrices
      mat_mult (mytempA, myB, myC, size);
			
      // Wait for receiving the "new" B matrix to arrive
      MPI_Wait (&req_recv, &status);
      MPI_Wait (&req_send, &status);

      memcpy (&myB[0], &mytempB[0], size * size * sizeof (double));
    }

  MPI_Barrier (MPI_COMM_WORLD);

	//Build the resultant matrix C
  MPI_Gatherv (&(myC[0]), size * size, MPI_DOUBLE, C, sendcounts, displs, blocktype, 0, comm_grid);

  MPI_Barrier (MPI_COMM_WORLD);
  if (grid_rank == 0)
    {
      timer = MPI_Wtime () - timer;
      all_time = MPI_Wtime() - start_time;
      printf ("%lf\n",timer);

      FILE * fp;
      fp = fopen("fox_time.txt","a");
      fprintf (fp, "%d, %d, %.8f, %.8f, %.8f, %.8f \n", N, p, timer, read_file_time, prep_time, all_time );
      fclose(fp);
    }
  mem_reset (&myA);
  mem_reset (&myB);
  mem_reset (&myC);
  mem_reset (&mytempA);
  mem_reset (&mytempB);
	///*
  if (grid_rank == 0)
    {
      //printf ("\nC=");
      //mat_print (C, N);
      mat_fprint (C, N, argv[2]);//comment this line if you dont want to write result to file
      mem_reset (&C);
    }
	//*/
  MPI_Type_free (&blocktype);
  MPI_Comm_free (&comm_row);
  MPI_Comm_free (&comm_col);
  MPI_Comm_free (&comm_grid);

  MPI_Finalize ();
  return 0;
}

// Dynamically allocate memory for an array
void *
mem_alloc (int n)
{
  void *ptr_mem = malloc (n);
  if (ptr_mem == NULL)
    {
      printf ("\nERROR! Not enough memory available.");
      MPI_Finalize ();
      exit (EXIT_FAILURE);
    }

  return ptr_mem;
}

// Free memory of dynamically allocated array
void
mem_reset (double **ptr)
{
  if (ptr)
    {
      free (*ptr);
      *ptr = NULL;
    }
}

// Fill matrix with a random double
void
mat_fill (double *A, int sz, FILE* fptr)
{
  for (int i = 0; i < sz * sz; i++)
    fscanf(fptr, "%lf", &A[i]);
}

// Multiply matrix blocks
void
mat_mult (double *A, double *B, double *C, int sz)
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

void
// Print matrix
mat_print (double *A, int sz)
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

void
// Print matrix to file
mat_fprint (double *A, int sz, char *outFileName)
{
  FILE *outFilePtr;
  outFilePtr = fopen(outFileName, "w");

  for (int i = 0; i < sz; i++)
    {
      for (int j = 0; j < sz; j++)
        fprintf (outFilePtr, "%.6f ", A[i * sz + j]);
      fprintf (outFilePtr,"\n");
    }
  fclose(outFilePtr);
}
