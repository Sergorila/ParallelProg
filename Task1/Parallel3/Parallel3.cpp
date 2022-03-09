#include <stdio.h>
int main(int argc, char* argv[])
{
	int count = 0;
#pragma omp parallel reduction (+: count)
	{
		count++;
		printf("Current value of the count: %d\n", count);
	}
	printf("Threads number: %d\n", count);
}