/*****************************************************
* Run by typing: "mpirun -np p ./qsort inputfile outputfile pivot_option"*
* p: Number of processors (square number)            *
*****************************************************/

#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(){
	int* a = NULL;
	a = (int*)malloc(2*sizeof(int));
	a = realloc(a, 0*sizeof(int));

	free(a);
	printf("%d\n",a[0]);
	printf("freed a\n");
	return 0;
}