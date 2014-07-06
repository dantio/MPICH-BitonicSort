#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

int N = 32; //N muss durch 8 Teilbar sein ansonsten muss man den algo auffüllen mit nullen

int *a;


void swap(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}


void bitonic(){

unsigned int i,j,k;
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<N; i++) {
		int ij=i^j;
		if ((ij)>i) {
			if ((i&k)==0 && a[i] > a[ij]) 
			{
				swap(i,ij);
			}
	      
			if ((i&k)!=0 && a[i] < a[ij])
			{
				swap(i,ij);
			}
		}
      }
    }
  }
	
}

void writeSorted(){
	int i;
	for (i=0;i < N;i++) {
		printf("%d ",a[i]);
	}
	printf("\n");
}

int main(int argc, char **argv) {
	
	a = malloc(N);
		
	int i;
	for (i = 0; i < N; i++) {
		a[i] = rand() % N;
     }
     
	bitonic();
	writeSorted();
	
}
