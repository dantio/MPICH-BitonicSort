/*
 ============================================================================
 Name        : ParallelBitonicSort.c
 Author      : Eric Wündisch, Daniil Tomilow
 Editversion BitonicSort from: https://www.cs.duke.edu/courses/fall08/cps196.1/Pthreads/bitonic.c
 ============================================================================
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <sys/time.h>

struct timeval startwtime, endwtime;
double seq_time;


int N;
int *a;

const int ASCENDING  = 1;
const int DESCENDING = 0;


void init(void);
void print(void);
void sort(void);
void test(void);
inline void swap(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void impBitonicSort(void);

int main(int argc, char **argv) {

  if (argc != 2) {
    printf("Usage: %s n\n  where n is problem size (power of two)\n", argv[0]);
    exit(1);
  }

  N = atoi(argv[1]);
  a = (int *) malloc(N * sizeof(int));

  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6+ endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n", seq_time);

  test();

  print();
}

void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "pass" : "fail");
}

void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N;
    //a[i] = (N - i);
  }
}

void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}

inline void swap(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}

void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j])) 
    swap(i,j);
}

void bitonicMerge(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++){
      compare(i, i+k, dir);
    }
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);
  }
}

void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
  }
}

void sort() {
  recBitonicSort(0, N, ASCENDING);
}


