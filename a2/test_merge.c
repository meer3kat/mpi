/*****************************************************
* Run by typing: "mpirun -np p ./qsort"    		 *
* p: Number of processors (square number)            *
*****************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


void local_merge(int size1, int size2, int* arr1, int* arr2, int* c){//allocate c before 
	// int size3=size1+size2;
	// int* c;
	// c = realloc(c, sizeof(int)*size3);
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


int main(){
	int len1 = 5;
	int len2 = 5;
	int* a1;
	int* a2;

	a1 = (int*)malloc(5*sizeof(int));
	a2 = (int*)malloc(5*sizeof(int));

	for(int i=0; i<5; i++){
		a1[i]=i+1;
		a2[i]=i*i-1;
		printf("%d, %d\n",i,a1[i]);
	}
	int* c;

	c=(int*)malloc(10*sizeof(int));
	local_merge(len1,len2,a1, a2,c );

	for(int i=0; i<10;i++){
		printf("%d, %d\n",i,c[i]);
	}
	
	return 0;
}



