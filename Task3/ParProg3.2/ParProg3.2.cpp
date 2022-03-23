#include <stdio.h>
#include <omp.h>
#include <locale.h>
int n;
#pragma omp threadprivate(n)
int main(int argc, char* argv[])
{
	n = 123;
#pragma omp parallel num_threads(4) copyin(n)
	{
		printf("Stream %d\n", omp_get_thread_num());
		printf("n: %d\n", n);
	}
}
