/*****************************************************
* Run by typing: "mpirun -np p ./qsort inputfile outputfile pivot_option"*
* p: Number of processors (square number)            *
*****************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void print_array(int* a, int length){
	printf("printing array below:\n");
	for(int i=0; i<length; i++){
		printf("%d. %d\n", i, a[i]);
	}
	printf("\nfinished\n");
}

int check_result(int *arr, int length)
{
    int i;
    for (i = 1; i < length; i++){
        if (arr[i - 1] > arr[i]){
            printf("error: a[%d] = %d, a[%d] = %d\n", i-1, arr[i-1], i, arr[i]);
            return -1;
        }
    }
    return 1;
    //printf("result correct\n");
}

void local_merge(int size1, int size2, int* arr1, int* arr2, int* c){//allocate c before 
	int i=0;
	int j=0;
	int k=0;
	if(arr1[size1-1]<arr2[size2-1]){
		while(j<size2){
			while(i<size1){
				if(arr1[i]<arr2[j]){
					c[k] = arr1[i];
					k++;
					i++;
				}
				else{
					c[k] = arr2[j];
					k++;
					j++;
				}
			}
			c[k]=arr2[j];
			k++;
			j++;
		}
	}
	else{
		while(i<size1){
			while(j<size2){
				if(arr1[i]<arr2[j]){
					c[k] = arr1[i];
					k++;
					i++;
				}
				else{
					c[k] = arr2[j];
					k++;
					j++;
				}
			}
			c[k]=arr1[i];
			k++;
			i++;
		}
	}

}


int partition(int* arr, int left, int right, int option){
	int i = left;
	// int j = right;
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

/* 
 * Copy d1 elements from index start with length steps into vector d2
 */ 
void copyArray(int *d1, int *d2, int start, int length) {
  int i;
  int j = start;
  for (i = 0; i < length; i++) {
    d2[i] = d1[j];
    j++;
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



void mpi_qsort(int* data, int len, MPI_Comm com, int option){
	MPI_Status status;
	int size, rank;
	MPI_Comm_size(com, &size);
	MPI_Comm_rank(com, &rank);
	MPI_Request req;
	int pivot;
	int* data_neighbour;
	int num_neighbour=0;
	int len_lo,len_hi;
	int* data_lo;
	int* data_hi;

	if(size == 1){
		// *last_length = len;
		// MPI_Send(data,len, MPI_INT, 0, 444, MPI_COMM_WORLD);
		MPI_Isend(data,len, MPI_INT, 0, 444, MPI_COMM_WORLD, &req);
		MPI_Request_free(&req);
		return;
	}

	if(option==0){
		if(rank == 0){pivot = data[len/2];} //set pivot to the middle of processoe 0. 
		MPI_Bcast(&pivot, 1, MPI_INT, 0, com);
	}
	else if(option==1){
		int processor_median = data[len/2];
		int *mean_median = NULL;
		if (rank == 0) {
  			mean_median = malloc(sizeof(int) * size);
		}
		MPI_Gather(&processor_median, 1, MPI_INT, mean_median, 1, MPI_INT, 0, com);
		if(rank ==0){
			// print_array(mean_median,size);
			quicksort(mean_median,0, size, option);
			pivot = mean_median[size/2];
			free(mean_median);
			mean_median = NULL;
		}
		MPI_Bcast(&pivot, 1, MPI_INT, 0, com);
	}
	else{
		int processor_median = data[len/2];
		// printf("processor median , %d\n", processor_median);
		int *mean_median;
		if (rank == 0) {
  			mean_median = malloc(sizeof(int) * size);
		}
		MPI_Gather(&processor_median, 1, MPI_INT, mean_median, 1, MPI_INT, 0, com);
		if(rank ==0){
			long int median_average = 0;
			for (int k=0; k<size; k++){
				// printf("pivot mean , %d\n", mean_median[k]);
				median_average = median_average + mean_median[k];

			}
			//print_array(mean_median,size);
			// printf("pivot mean , %ld\n", median_average);

			pivot = median_average/size;
			free(mean_median);
			mean_median = NULL;
			
		}
		MPI_Bcast(&pivot, 1, MPI_INT, 0, com);
		// printf("pivot, %ld\n", pivot);

	}


	int i = 0;
	while(i<len && data[i]<pivot){i++;}
	len_lo = i;
	len_hi = len-i;

	data_lo = (int*)malloc(len_lo*sizeof(int));
	data_hi = (int*)malloc(len_hi*sizeof(int));

	for(int j=0;j<i;j++) {data_lo[j]=data[j];} //write data to the left part low 
	for(int j=i;j<len;j++) {data_hi[j-i]=data[j];} //write data to the right part high
	
	//below to exchange data:
	int len_new;
	if(rank < size/2){
		MPI_Isend(data_hi,len_hi, MPI_INT, rank+size/2, rank, com, &req);
		//MPI_Send(data_hi,len_hi,MPI_INT,rank+size/2, rank, com);
		MPI_Probe(rank+size/2, rank+size/2, com, &status);
		MPI_Get_count(&status, MPI_INT, &num_neighbour);
		data_neighbour = (int*)malloc(num_neighbour*sizeof(int));
		MPI_Recv(data_neighbour, num_neighbour, MPI_INT, rank+size/2,rank+size/2,com, MPI_STATUS_IGNORE);
		
		// int* tmp;
		// tmp = (int*)malloc(len_lo*sizeof(int));
		// copyArray(data, tmp, 0, len_lo);

		// tmp 
		// int* tmp;
		// tmp = data;
		data = realloc(data, (len_lo+num_neighbour)*sizeof(int));
		local_merge(len_lo, num_neighbour, data_lo, data_neighbour, data);
		len_new = len_lo + num_neighbour;
		// print_array(data, len_new);
	}
	else{

		MPI_Probe(rank-size/2,rank-size/2,com, &status);
		MPI_Get_count(&status, MPI_INT, &num_neighbour);
		data_neighbour = (int*)malloc(num_neighbour*sizeof(int));
		MPI_Recv(data_neighbour, num_neighbour, MPI_INT,rank-size/2, rank-size/2, com, MPI_STATUS_IGNORE);
		MPI_Isend(data_lo, len_lo, MPI_INT, rank-size/2, rank, com, &req);	
		// MPI_Send(data_lo,len_lo,MPI_INT,rank-size/2, rank, com);
		// int* tmp;
		// tmp = data;
		data = realloc(data, (len_hi+num_neighbour)*sizeof(int));
		local_merge(len_hi, num_neighbour, data_hi, data_neighbour, data);
		len_new = len_hi + num_neighbour;
		// print_array(data, len_new);

	}

	MPI_Wait(&req, &status);
	free(data_lo);
	free(data_hi);
	free(data_neighbour);

	MPI_Comm sub;
	int color = rank/(size/2);
	MPI_Comm_split(com, color, rank, &sub);
	int n_size;
	MPI_Comm_size(sub, &n_size);
	mpi_qsort(data, len_new, sub, option);

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
	// printf("input file: %s\n",argv[0]);		
	// printf("pivot strategy %s\n", argv[3]);
	char* input_file = argv[1];
	char* output_file = argv[2];
	int option = atoi(argv[3]); //option number 
	// printf("%s\n",input_file);
	int* arr;   
	int n2;                          // create a pointer to the binary file data
	n2 = read_file(input_file, &arr);
	//**********************************************
	//read file and get n2, and chop it locally according to need. 

	// printf("rank: %d, n2 %d\n \n",rank, n2);
	// print_array(arr, n2);

	int chunk;             /* This many iterations will I do */
  	int istart, istop;  /* Variables for the local loop   */

	chunk  = n2/size;       /* Number of intervals per processor */
	istart = rank*chunk;         /* Calculate start and stop indices  */
	istop  = (rank+1)*chunk-1;     /* for the local loop                */

	if (rank == size-1 ) {
	istop = n2-1;           /* Make sure the last processor gets all the rest      */
	}
	int local_size;
	local_size = istop - istart + 1;
	int* local_arr = NULL;
	local_arr = (int*)malloc(local_size*sizeof(int));
	int local_index = 0;
	for(int i = istart; i<=istop; i++){
		local_arr[local_index] = arr[i];
		local_index++;
	}
	free(arr);
	// if(rank==1)print_array(arr,n2);



	//we read everything on every processor and assign work for them. 
	double t;
	if(rank == 0)
		t = MPI_Wtime ();
   	printf("I am before local quicksort\n");
	quicksort(local_arr,0,local_size-1,1); //local quick sort
	printf("finished  local quicksort, option: %d\n", option);
	// print_array(local_arr, local_size);
	//local sorted successfully. 
	MPI_Barrier(MPI_COMM_WORLD); 

	mpi_qsort(local_arr, local_size, MPI_COMM_WORLD,option);
	//and switch switch switch get mpi and ready to merge. 



	// MPI_Barrier(MPI_COMM_WORLD); 
	// if(rank==1)print_array(local_arr, local_size);
		// MPI_Barrier(MPI_COMM_WORLD); 


	int k=0;
	int num_get=0;
	int num_tmp;
	int collect_done = 0;
	int* sorted_array = NULL;

	int result=0;
	if(rank==0){
		// int* sorted_array;
		sorted_array = (int*)malloc(n2*sizeof(int));
		// int* len_final;
		// len_final = (int*)malloc(size*sizeof(int));
		while(k<size){
			MPI_Probe(k, 444, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &num_tmp);
			printf("num_tmp: %d from processor: %d", num_tmp, k);
			// MPI_Wait
			MPI_Recv(&sorted_array[num_get],num_tmp, MPI_INT, k, 444, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			num_get = num_get+ num_tmp;
			k++;
			// MPI_Barrier(MPI_COMM_WORLD); 
		}

	// 	// print_array(len_final,size);
		result = check_result(sorted_array,n2);
		collect_done = 1;
		MPI_Bcast(&collect_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD); 

	// 	// print_array(sorted_array,n2);
		// free(sorted_array);
	}
	else{
		MPI_Bcast(&collect_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
		//MPI_Barrier(MPI_COMM_WORLD); 

	}
	if(rank == 0) {
			t = MPI_Wtime () -t ;
			printf ("%d, %.8f, %d, %d, %d \n", n2, t, size, result, option);

			FILE * fp;
			fp = fopen ("A2outputb.txt","a");
			fprintf (fp, "%d, %.8f, %d, %d, %d \n", n2, t, size, result, option);
			fclose(fp);
	}

	//MPI_Bcast(sorted_array, n2, MPI_INT, 0, MPI_COMM_WORLD);

	if(rank==0){
		// print_array(sorted_array,n2);
		save_result(output_file, sorted_array, n2);
		free(sorted_array);	
		sorted_array = NULL;
	}
	printf("finished saving\n");
	MPI_Barrier(MPI_COMM_WORLD); 



	printf("local array: %d, rank: %d", local_arr[0],rank);

	if(local_arr != NULL) {
		printf("collected done from %d\n",rank);
		free(local_arr);
		printf("we freed local_arr from %d\n", rank);
		local_arr = NULL;
	}
	else{local_arr = NULL;}


	// if(rank==1){print_array(local_arr, local_size);}


		
	MPI_Finalize(); /* Shut down and clean up MPI */

	return 0;

}
