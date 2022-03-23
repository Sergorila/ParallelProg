#include <stdio.h>
#include <omp.h>
#include <locale.h>
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, ".ACP");
	int n;
#pragma omp parallel num_threads(3) private(n) //собственный счётчик для каждого потока
	{
		n = omp_get_thread_num();
		printf("Значение n (начало): %d\n", n);
#pragma omp single copyprivate(n) //теперь счётчик общий для всех потоков
		{
			n = 100;
		}
		printf("Значение n (конец): %d\n", n);
	}
}
