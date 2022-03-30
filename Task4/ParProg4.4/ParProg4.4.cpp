#include <stdio.h>
#include <omp.h>
#define MAX 8

int main() {
	int count = 0;
#pragma omp parallel num_threads(MAX)
	{
#pragma omp atomic //ATOMIC действует аналогично директиве CRITICAL, но относится только к идущему непосредственно за ней атомарному оператору 
		count++;
	}
	printf("Number of threads: %d\n", count);
}