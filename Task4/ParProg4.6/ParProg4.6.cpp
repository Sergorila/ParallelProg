#include <stdio.h>
#include <omp.h>
#include <locale.h>
#ifdef _WIN32
#include <windows.h>
void sleep(int Secnds)
{
	Sleep(100 * Secnds);
}
#endif
omp_lock_t lock;
int main(int argc, char* argv[])
{
	setlocale(LC_ALL, ".ACP");
	int n;
	omp_init_lock(&lock);
	printf("Блокировка инициализированна\n");
#pragma omp parallel num_threads(4) private (n)
	{
		n = omp_get_thread_num();
		while (omp_test_lock(&lock) == 0)
		{
			printf("Секция закрыта, поток %d\n", n);
			sleep(2);
		}
		printf("Начало закрытой секции, поток %d\n", n);
		sleep(2 + 3 * ((n + 1) % 2));
		printf("Конец закрытой секции, поток %d\n", n);
		omp_unset_lock(&lock);
	}
	omp_destroy_lock(&lock);
	printf("Блокировка не инициализированна\n");
}
