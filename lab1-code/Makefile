############################
# Makefile for lab 1
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -O2
LIBS = -lm

BINS = hello deadlock exchange pingpong onetoall pi buffon-laplace

all: $(BINS)

buffon-laplace: buffon-laplace.c buffon-laplace.h
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

