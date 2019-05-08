/*****************************************************
* Run by typing: "mpirun -np p ./qsort"    		 *
* p: Number of processors (square number)            *
*****************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void exchange(int* a, int* b){
    int tmp = *a;
    *a = *b;
    *b = tmp;
}

void check_result(int *arr, int length)
{
    int i;
    for (i = 1; i < length; i++){
        if (arr[i - 1] > arr[i]){
            printf("error: a[%d] = %d, a[%d] = %d\n", i-1, a[i-1], i, a[i]);
        }
    }
}

int partition(int* arr, int left, int right, int option){
	int i = left;
	int j = right;
	int tmp;
	int pivot;
	if(option == 0){
		pivot = arr[right];
	}
	else if(option == 1){
		pivot = arr[(left+right)/2];
		arr[(left+right)/2] = arr[right];
		arr[right] = pivot;
	}
	else{
		pivot = arr[right];
	}	
	//printf("i: %d \t, j: %d \t, pivot index: %d \t, pivot: %d \t",i,j,(left+right)/2,pivot);
	i = left - 1;
	for(int k = left; k < right; k++){
		if (arr[k] <= pivot){
			i++;
			tmp = arr[i];
			arr[i] = arr[k];
			arr[k] = tmp;
		}
	}
	tmp = arr[i+1];
	arr[i+1]=arr[right];
	arr[right] = tmp;

	return i+1;
}

void quicksort(int* arr, int left, int right, int option){
	if (left < right){
		int pindex = partition(arr, left, right, option);
		//printf("pindex %d,left %d, right %d\n",pindex, left, right);
		if(left<pindex) { quicksort(arr,left, pindex-1,option);}
		if(pindex+1<right) { quicksort(arr, pindex+1, right,option);}
	}
}

void mpi_qsort(int* begin, int loc_size, MPI_Comm comm){
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
		int pivot_index=0;
		if(rank==0){
			pivot_in
		}	
	}
}

// void exchange(int p1, int p2, int size1, int size2, int* arr1, int* arr2){
// 	int pivot1, pivot2;
// 	pivot1 = size1/2;
// 	pivot2 = size2/2;
// 	MPI_Send(&arr1,)
// }


void local_merge(int size1, int size2, int* arr1, int* arr2, int* c){//allocate c before 
	int size3;
	// int* c;
	// c = (int*)malloc(sizeof(int)*size3);
	int i=0;
	int j=0;
	int k=0;
	if(arr1[size1-1]<arr2[size2-1]){
		while(j<size2){
			while(i<size1){
				if(arr1[i]<arr2[j]){
					c[k] = arr1[i];
					i++;
				}
				else{
					c[k] = arr2[j];
					j++;
				}
			}
			c[k]=arr2[j];
			j++;
		}
	}
	else{
		while(i<size1){
			while(j<size2){
				if(arr1[i]<arr2[j]){
					c[k] = arr1[i];
					i++;
				}
				else{
					c[k] = arr2[j];
					j++;
				}
			}
			c[k]=arr1[i];
			i++;
		}
	}

}


int read_file(char *name, int** pp){
	//return the number of number read in the file and pointer *pp will point to the first
    FILE* f;
    f = fopen(name, "r");

    if(f){
        //printf("file opened! \n");
        fseek(f, 0, SEEK_END);
        fseek(f, 0, SEEK_SET);
        int *p = NULL;
        int n;
        fscanf(f,"%d ",&n);
        p = (int*)malloc(n*sizeof(int));
        for(int i=0;i<n;i++){
        	fscanf(f,"%d ", &p[i]);
        }
        *pp = &p[0];
        fclose(f);
        return n;
    }
    else
    	return 0;
}

void save_result(char* name, int* arr, int n){
	FILE* f;
	f = fopen(name,"w");
	for(int i=0; i<n; i++){
		fprintf(f,"%d\n",arr[i]);
	}
	fclose(f);
}


int main(int argc, char *argv[]){	
	// set up
	MPI_Init(&argc, &argv); //initialize 
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //get my number
	MPI_Comm_size(MPI_COMM_WORLD, &size); //get number of processors
	MPI_Status status;
	// Program arguments
	if(argc != 4){
		printf("please enter 3 input; to run ./qsort inputfile outputfile methods\n.");
		return 1;
	}
	printf("yes.\n");
	printf("input file: %s\n",argv[0]);		
	printf("pivot strategy %s\n",argv[2]);
	char* input_file = argv[1];
	char* output_file = argv[2];
	int option = atoi(argv[3]); //option number 
	printf("%s\n",input_file);
	int* arr;   
	int n2;                          // create a pointer to the binary file data
	n2 = read_file(input_file,&arr);
	printf("rank: %d, n2 %d\n \n",rank, n2);

	for(int i=0;i<n2;i++){
		printf("%d. %d \n",i, arr[i]);
	}
	int chunk;             /* This many iterations will I do */
  	int i, j, istart, istop;  /* Variables for the local loop   */

	chunk  = n2/size;       /* Number of intervals per processor */
	istart = rank*chunk;         /* Calculate start and stop indices  */
	istop  = (rank+1)*chunk-1;     /* for the local loop                */

	if (rank == size-1 ) {
	istop = n2-1;           /* Make sure the last processor gets all the rest      */
	}
	int local_size;
	local_size = istop - istart + 1;
	int* local_arr;
	local_arr = (int*)malloc(local_size*sizeof(int));
	int local_index = 0;
	for(int i = istart; i<=istop; i++){
		local_arr[local_index] = arr[i];
		local_index++;
	}


	free(arr);           
	quicksort(local_arr,0,local_size-1,option);                


	printf("result\n");


	for(int i=0;i<local_size;i++){
		printf("%d. %d \n",i, local_arr[i]);
	}

	//save_result(output_file, arr, n2);
	printf("finished saving\n");
	free(local_arr);	
	MPI_Finalize(); /* Shut down and clean up MPI */	
	return 0;

}
