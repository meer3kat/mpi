
void mpi_qsort(int* begin, int loc_size, MPI_Comm comm, int* last_length){
	int comm_size;
	int rank;
	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &rank);
	int big_size = comm_size+2*loc_size;
	int* loc_array = (int*)malloc(sizeof(int) * big_size);

	printf("current size is %d comm size is %d \n", loc_size, comm_size);
	int* initial_sizes = (int*)malloc(sizeof(int) * comm_size);
	MPI_Allgather(&loc_array, 1, MPI_Int, &initial_sizes[0], 1, MPI_Int,comm);
	int i=0;
	for(i=0;i<local_size;i++){
		loc_array[i] = begin[i];
	}
	recursive_sort(loc_array,loc_size,comm);

	int* sorted_sizes = (int*)malloc(sizeof(int) * comm_size);
	MPI_Allgather(&local_size, 1, MPI_Int, &sorted_sizes[0], 1, MPI_Int, comm);

	int* sendcounts = (int*)malloc(sizeof(int) * comm_size);
    int* senddispls = (int*)malloc(sizeof(int) * comm_size);
    int* recvcounts = (int*)malloc(sizeof(int) * comm_size);
    int* recvdispls = (int*)malloc(sizeof(int) * comm_size); 


    // calculate sendcounts, senddispls, recvcounts, recvdispls to send it to MPI_Alltoallv to perform merging
    Set_Sending_Size(comm_size, rank, sorted_sizes, initial_sizes, sendcounts, senddispls ); 
	Set_Receiving_Size(comm_size, rank, sorted_sizes, initial_sizes, recvcounts, recvdispls); 
    // Sends data from all to all processes
    MPI_Alltoallv(internal_array, sendcounts, senddispls, MPI_INT, begin, recvcounts, recvdispls, MPI_INT, comm); 


	// free all allocated temporary arrays
    MPI_Barrier(comm); 
    free(internal_array); 
    free(initial_sizes); 
    free(sorted_sizes); 
    free(sendcounts); 
	free(senddispls);
    free(recvcounts); 
    free(recvdispls); 
	return ;
}

void recursive_sort(int* begin, int local_size, MPI_Comm comm){
	int comm_size;
	int rank;
	MPI_Comm_size(comm, &comm_size);
	MPI_Comm_rank(comm, &rank);
	if(comm_size>1){
		int tmp = 0;
		int i =0;

		int root_size = 0;
		MPI_Allreduce(&local_size, &root_size, 1, MPI_INT, MPI_SUM, comm);
		int pivot_value=0;
		if(rank==0){
			pivot_value = begin[local_size/2];
		}
		MPI_Bcast(&pivot_value,1,MPI_INT, 0, MPI_COMM_WORLD); //brocast the pivot to everyone 

	
	}
}