#include <stdio.h>
#include <omp.h>
#include <locale.h>
#define NNN 10
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, ".ACP");
	int i, n;
#pragma omp parallel private (i, n)
	{
		n = omp_get_thread_num();
#pragma omp for ordered
		for (i = 0; i < NNN; i++)
		{
#pragma omp ordered
			{
				printf("Поток %d, итерация %d\n", n, i);
			}
		}
	}
}