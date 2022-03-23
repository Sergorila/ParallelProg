#include <stdio.h>
#include <omp.h>
#include <locale.h>
#define size 35
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, ".ACP");
	int A[size], B[size], C[size], i, n;
	for (i = 0; i < size; i++) { A[i] = i; B[i] = 2 * i; C[i] = 0; }
#pragma omp parallel num_threads(4) shared(A, B, C) private(i, n)
	{
		n = omp_get_thread_num();
#pragma omp for
		for (i = 0; i < size; i++)
		{
			C[i] = A[i] + B[i];
			printf("Поток %d сложил элементы с номером %d\n",
				n, i);
			printf("%d + %d = %d\n", A[i], B[i], C[i]);
		}
	}
}