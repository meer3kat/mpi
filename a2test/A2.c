//Assignment2: Lomuto partition scheme Quicksort using MPI

#include <mpi.h>
#include <stdio.h>
#include<stdlib.h>
#include <sys/time.h>
#include<time.h>
#include <string.h> 
#define print_elements 1


void quicksort(double* A, int lo, int hi);
double partition(double* A, int lo, int hi);
void pquick(double* data,int len, MPI_Comm com);
void merge(double* a, int m, double* b, int n, double* sorted);
void print_array(double *A, int len);


int read_file(char *name, double** pp){
	//return the number of number read in the file and pointer *pp will point to the first
    FILE* f;
    f = fopen(name, "r");

    if(f){
        //printf("file opened! \n");
        fseek(f, 0, SEEK_END);
        fseek(f, 0, SEEK_SET);
        double *p = NULL;
        int n;
        fscanf(f,"%d ",&n);
        p = (double*)malloc(n*sizeof(double));
        for(int i=0;i<n;i++){
        	fscanf(f,"%lf ", &p[i]);
        }
        *pp = &p[0];
        fclose(f);
        return n;
    }
    else
    	return 0;
}


int main(int argc, char *argv[]){
  /*check input*/
  if (argc != 2) {
    printf("Input error.Please input the number of Elements.");
    exit(EXIT_FAILURE);
  }
  int num_element;
  // char* input_file = argv[1];


	// int n2;                          // create a pointer to the binary file data
	// n2 = read_file(input_file,&arr);
	// printf("rank: %d, n2 %d\n \n",rank, n2);

  // num_element= read_file(input_file, &)
  

  srand48(time(NULL));
  int rank,size,rc;
	double start,end;//record time

	/*Initialize MPI*/
  int check_mpi=MPI_Init(&argc, &argv);
  if(check_mpi != MPI_SUCCESS) {
    printf("Error in starting MPI program. \n");
    MPI_Abort(MPI_COMM_WORLD, check_mpi );
  }
	start = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &size); /* Get the number of processors */
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number*/

  double *A=NULL;
  double *sub_A=NULL;

/*Generate data on P0*/
  if(rank==0){

	char* input_file = argv[1];
	printf("%s\n",input_file);
	double* A;   
	int n2;                          // create a pointer to the binary file data
	num_element = read_file(input_file,&A);
	printf("rank: %d, n2 %d\n \n",rank, num_element);

	
#if print_elements
    printf("Elements:");
    print_array(A,num_element);
#endif

  }

/*Scatter data and Sort locally*/
  int sub_num=num_element/size;
  sub_A=(double *)malloc(sub_num*sizeof (double));
  MPI_Scatter(A, sub_num, MPI_DOUBLE,sub_A , sub_num,MPI_DOUBLE, 0,MPI_COMM_WORLD);
  quicksort(sub_A,0,sub_num-1);
  MPI_Barrier(MPI_COMM_WORLD);
#if print_elements
  printf("rank %d locally sorted: ",rank);
  print_array(sub_A, sub_num);
#endif

/*Call pquick*/
  pquick(sub_A,sub_num, MPI_COMM_WORLD);

/*P0 probe and receive data from all other*/
  MPI_Status status;
  int k=0;
  int num_get=0;
  int num_tmp;
  if(rank==0 ){
	  while(k<size){
  	  MPI_Probe(k, 444,MPI_COMM_WORLD, &status);
		  MPI_Get_count(&status, MPI_DOUBLE, &num_tmp);
		  MPI_Recv(&A[num_get],num_tmp,MPI_DOUBLE,k,444
,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		  num_get=num_get+num_tmp;
		  k++;
	  }
#if print_elements
	 printf("Sorted Elements:");
   print_array(A,num_element);
#endif
	 end = MPI_Wtime();
	 printf("proc_num:%d, num_element:%d, time:%f\n", size, num_element, end-start);	
  }

/*check result*/
  if(rank==0){
    for(int i=0;i<num_element-1;i++){
      if(A[i]>A[i+1]){
        printf("result error!\n");
        free(A);
        return 1;
       }
    }
    printf("The array has been sorted!\n");
		free(A);
  }
  //if(sub_A!=NULL) free(sub_A);
  MPI_Finalize();
  return 0;
}

void quicksort(double* A, int lo, int hi){
	double pivot;
	if(lo<hi){
		pivot=partition(A,lo,hi);
		quicksort(A,lo,pivot-1);
		quicksort(A,pivot+1,hi);
	}
}

double partition(double* A,int lo, int hi){

	double temp;
	int i=lo-1;
	for(int j=lo;j<hi;j++){
		if(A[j]<A[hi]){
			i++;
			temp=A[i];A[i]=A[j];A[j]=temp;
		}
	}
	
	temp=A[i+1];A[i+1]=A[hi];A[hi]=temp;
	return i+1;
}

void pquick(double* data,int len, MPI_Comm com){

	MPI_Status status;
	int size,rank;
	MPI_Comm_size(com, &size);
	MPI_Comm_rank(com, &rank);
	MPI_Request req;
	double pivot;
	double * data_neighbour;
	int num_neighbour=0;
	int len_lo,len_hi;
	double* data_lo;
	double* data_hi;

	if(size==1){
		MPI_Isend(data,len, MPI_DOUBLE,0,444,MPI_COMM_WORLD,&req);
		MPI_Request_free(&req);
		 return;
		}

	if(rank==0){
		pivot=data[len/2];
		
	}
	MPI_Bcast(&pivot,1,MPI_DOUBLE,0,com);

/*Split data*/
	int i=0;
	while(i<len && data[i]<pivot) i++;
	len_lo=i;
	len_hi=len-i;

	data_lo=(double *) malloc (i*sizeof (double));
	data_hi=(double *) malloc ((len-i)*sizeof (double));
	for(int j=0;j<i;j++){
		data_lo[j]=data[j];
	}
	for(int j=i;j<len;j++){
		data_hi[j-i]=data[j];
	}

/*Exchange data pairwise and merge*/
  int len_new;
  if(rank<size/2){
		MPI_Isend(data_hi,len_hi, MPI_DOUBLE,rank+size/2,rank,com,&req);
		MPI_Probe(rank+size/2,rank+size/2,com, &status);
    MPI_Get_count(&status, MPI_DOUBLE, &num_neighbour);

    data_neighbour=(double *) calloc (num_neighbour,sizeof (double));
    MPI_Recv(data_neighbour,num_neighbour,MPI_DOUBLE,rank+size/2,rank+size/2
,com, MPI_STATUS_IGNORE);
	//merge
		double* tmp;
		tmp=data;
		data=realloc(tmp,(len_lo+num_neighbour)*sizeof(double));
		merge(data_lo,len_lo, data_neighbour,num_neighbour,data);
		len_new=len_lo+num_neighbour;
#if print_elements
		printf("rank:%d ",rank);
		print_array(data,len_lo+num_neighbour);
#endif

	}else{

	
		MPI_Isend(data_lo,len_lo, MPI_DOUBLE,rank-size/2,rank,com,&req);
    MPI_Probe(rank-size/2,rank-size/2,com, &status);
		MPI_Get_count(&status, MPI_DOUBLE, &num_neighbour);
	  data_neighbour=(double *) malloc (num_neighbour*sizeof (double));
    MPI_Recv(data_neighbour,num_neighbour,MPI_DOUBLE,rank-size/2,rank-size/2,com,&status);
		//merge
		double* tmp;
		tmp=data;
		data=realloc(tmp,(len_hi+num_neighbour)*sizeof(double));
		merge(data_hi,len_hi, data_neighbour,num_neighbour,data);
			len_new=len_hi+num_neighbour;
#if print_elements
		printf("rank:%d ",rank);
		print_array(data,len_hi+num_neighbour);
#endif
	
	}

	MPI_Wait(&req,&status);
	free(data_lo);
	free(data_hi);
	free(data_neighbour);

	MPI_Comm sub;
	MPI_Comm_split(com,rank<size/2,0,&sub);
	int n_size;
	MPI_Comm_size(sub, &n_size);   
	pquick(data,len_new,sub);	
}


void merge(double* a, int m, double* b, int n, double* sorted) {

	
  int i, j, k;
 
  j = k = 0;
 
  for (i = 0; i < m + n;) {
    if (j < m && k < n) {
      if (a[j] < b[k]) {
        sorted[i] = a[j];
        j++;
      }
      else {
        sorted[i] = b[k];
        k++;
      }
			i++;
    }
    else if (j == m) {
      while (i < m + n) {
        sorted[i] = b[k];
        k++;
        i++;
      }
    }
    else {
      while(i < m + n) {
        sorted[i] = a[j];
        j++;
        i++;
      }
    }
  }

}
void print_array(double *A, int len){
	int i;
	for(i=0;i<len;i++){
			printf("%lf ",A[i]);
		}
	printf("\n");
}



		  
