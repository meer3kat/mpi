/*****************************************************
* Run by typing: "mpirun -np p ./qsort"    		 *
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

void check_result(int *arr, int length)
{
    int i;
    for (i = 1; i < length; i++){
        if (arr[i - 1] > arr[i]){
            printf("error: a[%d] = %d, a[%d] = %d\n", i-1, arr[i-1], i, arr[i]);
            return;
        }
    }
    printf("result correct\n");
}



void local_merge(int size1, int size2, int* arr1, int* arr2, int* c){//allocate c before 
	int size3=size1+size2;
	// int* c;
	c = realloc(c, sizeof(int)*size3);
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
	// last_length = local_size;    
	quicksort(local_arr,0,local_size-1,option);

	print_array(local_arr, local_size);


	// mpi_qsort(local_arr, local_size, MPI_COMM_WORLD, last_length);  
	// local_size = *last_length;
	// MPI_Gather(&local_size, 1, MPI_INT, receive_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// if(rank == 0){
	// 	int index = 0;
	// 	receive_count[0]= index;
	// 	for(int i=1;i<size; i++){
	// 		index = index+receive_count[i-1];
	// 		receive_displacement[i] = index;
	// 	}
	// } 

	// MPI_Gatherv(local_arr,local_size, MPI_INT, data_sorted, receive_count,receive_displacement,MPI_INT,0, MPI_COMM_WORLD);

	// // printf("print local array\n");
	// // for(int i = 0; i< local_size; i++){
	// // 	printf("%d. %d \n",i, local_arr[i]);
	// // }

	// if(rank == 0){
	// 	printf("print result below\n");
	// 	for(int i=0;i<n2;i++){
	// 		printf("%d. %d \n",i, arr[i]);
	// 	}
	// 	printf("finished printing\n");
	// }

	// if(rank == size-1){
	// 	printf("checking result");
	// 	check_result(arr, n2);

	// }


	// printf("result\n");


	// for(int i=0;i<local_size;i++){
	// 	printf("%d. %d \n",i, local_arr[i]);
	// }

	//save_result(output_file, arr, n2);
	printf("finished saving\n");
	free(local_arr);
	free(arr);	
	MPI_Finalize(); /* Shut down and clean up MPI */	
	return 0;

}
