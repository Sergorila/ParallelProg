#include <iostream>
#include <omp.h>
using namespace std;
int main()
{
	double start_time = omp_get_wtime();
	double end_time = omp_get_wtime();
	double tick = omp_get_wtick();
	cout << end_time - start_time << endl;
	cout << "Accuracy timer" << tick << endl;
}
