#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "buffon-laplace.h"

/**
 * This program estimates the value of PI using a Buffon-Laplace simulation.
 * @author Malin Kallen 2019
 */
int main(int argc, char **argv) {
	// Set up
	MPI_Init(&argc, &argv);
	int rank, num_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	// Program arguments
	if (5  != argc) {
		if (0 == rank) {
			printf("Usage: ./buffon-laplace A B L N, where\n");
			printf(" A = horizontal spacing of grid lines\n");
			printf(" B = vertical spacing of grid lines\n");
			printf(" L = length of needle\n");
			printf(" N = total number of times the needle is dropped\n");
		}
		return 1;
	}

	const double A = atof(argv[1]);
	const double B  = atof(argv[2]);
	const double L = atof(argv[3]);
	const int globalN = atoi(argv[4]);
	const int localN = globalN/num_proc + (rank < globalN%num_proc);

	// Simulate needle throws and count number of crossings
	int localC = run_simulation(A, B, L, localN);

	// Approximate pi
	int globalC;
	MPI_Reduce(&localC, &globalC, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (0 == rank) {
		double pi_approx = globalN * (2*L*(A+B)-L*L) / (globalC*A*B);
		printf("N: %d, C: %d, PI: %f\n", globalN, globalC, pi_approx);
	}

	MPI_Finalize();
	return 0;
}

int run_simulation(double A, double B, double L, int N) {
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	srand(rank);	// Not random, but at least different in different processes!
	int C = 0;
	point_t eye, point;
	for (int i=0; i<N; i++) {
		throw_needle(A, B, L, &eye, &point);
		if (crossing(point, A, B)) {
			C++;
		}
	}
	return C;
}

void throw_needle(double A, double B, double L, point_t *eye, point_t *point) {
	eye->x = A * rand()/RAND_MAX;
	eye->y = B * rand()/RAND_MAX;
	double alfa = 2 * PI * rand()/RAND_MAX;	// Angle of needle
	point->x = eye->x + L * cos(alfa);
	point->y = eye->y + L * sin(alfa);
}

int crossing(point_t point, double A, double B) {
	return point.x<0 || point.x>A || point.y<0 || point.y>B;
}

