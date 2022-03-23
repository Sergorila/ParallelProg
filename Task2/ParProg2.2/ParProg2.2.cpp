#include <stdio.h>
#include <locale.h>
#include <omp.h>
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, ".ACP");
#pragma omp parallel num_threads(3)
	{
		printf("Печать сообщения 1\n");
#pragma omp single
		{
			printf("Выполняем в одном потоке\n");
		}
		printf("Печать сообщения 2\n");
	}
}