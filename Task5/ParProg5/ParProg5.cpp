#include <iostream>
#include <cmath>
#include <time.h>
#include <omp.h>
#define N 1000
using namespace std;
int main()
{
	int i, j, k;
	double h = 1.0 / N;
	double A = 0, B = 0;
	double Tms;
	Tms = clock();
	// основной вычислительный блок

#pragma omp parallel for num_threads(4) reduction(+: A) shared(i, j, k) 
	for (i = 0; i < N - 1; i++) {
		for (j = 0; j < N - 1; j++) {
			for (k = 0; k < N - 1; k++) {
				A += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
			}
		}
	}
	A *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time = " << Tms << " sec" << endl;
	cout << "A = " << A << endl;
	A = 0;


#pragma omp parallel for num_threads(4) reduction(+: A) schedule (static) shared(i, j, k)
		for (i = 0; i < N - 1; i++){
			for (j = 0; j < N - 1; j++) {
				for (k = 0; k < N - 1; k++) {
					A += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
				}
			}
		}
	A *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time static = " << Tms << " sec" << endl;
	cout << "A = " << A << endl;
	A = 0;

#pragma omp parallel for num_threads(4) reduction(+: A) schedule (dynamic) shared(i, j, k)
	for (i = 0; i < N - 1; i++) {
		for (j = 0; j < N - 1; j++) {
			for (k = 0; k < N - 1; k++) {
				A += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
			}
		}
	}
	A *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time dynamic = " << Tms << " sec" << endl;
	cout << "A = " << A << endl;
	A = 0;

#pragma omp parallel for num_threads(3) collapse(3) reduction(+: B)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N - 1; j++) {
			for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
				B += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
		}
	}
	B *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time = " << Tms << " sec" << endl;
	cout << "B = " << B << endl;
	B = 0;

#pragma omp parallel for num_threads(3) collapse(3) reduction(+: B) schedule (static)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N - 1; j++) {
			for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
				B += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
		}
	}
	B *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time static = " << Tms << " sec" << endl;
	cout << "B = " << B << endl;
	B = 0;

#pragma omp parallel for num_threads(3) collapse(3) reduction(+: B) schedule (dynamic)
	for (i = 0; i < N; i++) {
		for (j = 0; j < N - 1; j++) {
			for (k = 0; k < int(sqrt(pow((N - i), 2) - pow(j, 2))); k++)
				B += cos(h * (i + 1.0 / 2.0)) / (1 + exp(-(pow(h * (i + 1.0 / 2.0), 2) + pow(h * (j + 1.0 / 2.0), 2) + pow(h * (k + 1.0 / 2.0), 2))));
		}
	}
	B *= pow(h, 3);
	Tms = (clock() - Tms) / CLOCKS_PER_SEC;
	cout << "Time dynamic = " << Tms << " sec" << endl;
	cout << "B = " << B << endl;
	B = 0;
}