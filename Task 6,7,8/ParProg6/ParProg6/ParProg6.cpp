#include <mpi.h>
#include <iostream>
#include <stdlib.h>
#include <cmath>

using namespace std;

void task_6_2();
void task_6_4();
void task_6_6();
void task_7_2();
void task_7_4();


bool Is_Finished(bool* Prizn, int N)
{
	bool Res = true;
	for (int k = 0; k < N; k++)
		Res = (Prizn[k]) && Res;
	return Res;
}

double Master(double* A, int N)
{
	double S = 0;
	for (int k = 0; k < N; k++)
		S += A[k];
	return S;
}

void Slave(double* A, int N, int Rank)
{
	static int Count = 0;
	double* Tmp1 = new double[N];
	double* Tmp2 = new double[N];
	Tmp1[0] = cos(1.0 + Count + Rank);
	Tmp2[0] = cos(1.0 + Count - Rank);
	for (int k = 1; k < N; k++)
	{
		Tmp1[k] = sin(Tmp1[k - 1] + Tmp2[k - 1]);
		Tmp2[k] = sin(Tmp1[k - 1] - Tmp2[k - 1]);

	}
	for (int k = 0; k < N; k++) 
	{
		A[k] = 0;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				A[k] += Tmp1[i] * Tmp2[j] / (1 + k + (i - j) * (i - j));
	}
	Count++;
	delete[] Tmp2;
	delete[] Tmp1;
}

void task_7_6();
void task_8_2();
void task_8_4();
void task_8_6();

int main(int argc, char** argv)
{
	task_8_6();
	return 0;
}

void task_6_2()
{
	int const NTimes = 100;
	char Proc_Name[MPI_MAX_PROCESSOR_NAME + 1];
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int NLen;
	MPI_Get_processor_name(Proc_Name, &NLen);
	double tick = MPI_Wtick();
	double time_start = MPI_Wtime();
	double time_finish;
	for (int i = 0; i <= NTimes; i++)
	{
		time_finish = MPI_Wtime();
	}
	cout << "Processor" << rank << ' ' << Proc_Name << endl;
	cout << "timer's tick = " << tick << " time = " <<
		(time_finish - time_start) / NTimes << endl;
	MPI_Finalize();
}

void task_6_4()
{
	double a = 0;
	double b = 0;
	MPI_Init(NULL, NULL);
	int Size;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	int Rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	if (Rank == 0)
	{
		a = 7.7;
		MPI_Send(&a, 1, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);
		MPI_Status Status;
		MPI_Recv(&b, 1, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD, &Status);
	}
	if (Rank == 1)
	{
		b = 3.3;
		MPI_Status Status;
		MPI_Recv(&a, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &Status);
		MPI_Send(&b, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
	}
	cout << "Process " << Rank << " a=" << a << " b=" << b << endl;
	MPI_Finalize();
}

void task_6_6()
{
	MPI_Init(NULL, NULL);
	int Size;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	if (Size >= 2)
	{
		int Rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
		if (Rank == 0)
		{
			int N = 3 + (rand() % 10);
			double* A = new double[N];
			for (int k = 0; k < N; k++)
			{
				A[k] = k + 0.7;
			}
			cout << "Process 0 sends the data: " << endl;
			for (int k = 0; k < N; k++)
			{
				cout << A[k] << endl;
			}
			MPI_Send(A, N, MPI_DOUBLE, 1, 5, MPI_COMM_WORLD);
			delete[] A;
		}
		if (Rank == 1)
		{
			MPI_Status St;
			MPI_Probe(0, 5, MPI_COMM_WORLD, &St);
			int N;
			MPI_Get_count(&St, MPI_DOUBLE, &N);
			double* A = new double[N];
			MPI_Recv(A, N, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &St);
			cout << "Process 1 has accepted the data:" << endl;
			for (int k = 0; k < N; k++)
			{
				cout << A[k] << endl;
			}
			delete[] A;
		}
		else
			cout << "Performance of the program not less than in 2 processes is required" << endl;
		MPI_Finalize();
	}
}

void task_7_2()
{
	double a = 0;
	double b = 0;
	double c = 0;
	int const BufSize = sizeof(double) + MPI_BSEND_OVERHEAD;
	MPI_Init(NULL, NULL);
	int Size;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	int Rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	int Posl = (Rank + 1) % Size; //Номер последующего процесса
	int Pred = Rank ? Rank - 1 : Size - 1; //Номер предшествующего процесса
	a = Rank + 0.7;
	//void *Addr=NULL; 
	void *Buf = malloc(BufSize);
	MPI_Buffer_attach(Buf, BufSize);
	MPI_Bsend(&a, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD);
	MPI_Status St;
	MPI_Recv(&b, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &St);
	int Tmp; 
	MPI_Buffer_detach(Buf, &Tmp);
	MPI_Buffer_attach(Buf, BufSize);
	MPI_Bsend(&a, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD);
	MPI_Recv(&c, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &St);
	MPI_Buffer_detach(Buf, &Tmp);
	if (!Buf)
		free(Buf);
	cout << "Process " << Rank << " a=" << a << " b=" << b << " c=" << c << endl;
	MPI_Finalize();
}

void task_7_4()
{
	double a = 0;
	double b = 0;
	double c = 0;
	MPI_Init(NULL, NULL);
	int Size;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	int Rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	int Posl = (Rank + 1) % Size; //Номер последующего процесса
	int Pred = Rank ? Rank - 1 : Size - 1; //Номер предшествующего процесса
	a = Rank+0.7;
	MPI_Status St;
	MPI_Sendrecv(&a, 1, MPI_DOUBLE, Posl, 5, 
		&b, 1, MPI_DOUBLE, Pred, 5, MPI_COMM_WORLD, &St);
	MPI_Sendrecv(&a, 1, MPI_DOUBLE, Pred, 5,
		&c, 1, MPI_DOUBLE, Posl, 5, MPI_COMM_WORLD, &St);
	cout << "Process " << Rank << " a=" << a << " b=" << b << " c=" << c << endl;
	MPI_Finalize();
}
void task_7_6()
{
	int NN = 1000, NNN = 3;
	MPI_Init(NULL, NULL);
	int Size, Rank;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);

	if (Size < 2)
		std::cout << "Minimum two processes required." << std::endl;
	else
	{
		if (Rank == 0)
		{
			bool* Prizn = new bool[Size - 1];
			double* A = new double[(Size - 1) * NN];
			double* R = new double[Size - 1];

			for (int k = 0; k < Size - 1; ++k)
				R[k] = 0;

			MPI_Request* Request = new MPI_Request[Size - 1];
			MPI_Status* Status = new MPI_Status[Size - 1];

			int* Indx = new int[Size - 1];
			double Tms = MPI_Wtime();

			for (int m = 0; m < NNN; ++m)
			{
				for (int k = 1; k < Size; ++k)
				{
					Prizn[k - 1] = false;
					MPI_Irecv(&A[NN * (k - 1)], NN, MPI_DOUBLE, k, 5, MPI_COMM_WORLD, &Request[k - 1]);
				}

				while (!Is_Finished(Prizn, Size - 1))
				{
					int Num;
					MPI_Waitsome(Size - 1, Request, &Num, Indx, Status);

					for (int k = 0; k < Num; ++k)
					{
						Prizn[Indx[k]] = true;
						R[Indx[k]] += Master(&A[NN * Indx[k]], NN);
					}
				}
			}

			Tms = MPI_Wtime() - Tms;

			for (int k = 0; k < Size - 1; ++k)
				std::cout << "R[" << k << "]=" << R[k] << std::endl;

			std::cout << "Time: " << Tms << " s" << std::endl;

			delete[] Indx;
			delete[] Status;
			delete[] Request;
			delete[] R;
			delete[] A;
			delete[] Prizn;
		}
		else
		{
			double* A = new double[NN];

			for (int m = 0; m < NNN; ++m)
			{
				Slave(A, NN, Rank);
				MPI_Send(A, NN, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
			}

			delete[] A;
		}
	}
	MPI_Finalize();
}

void task_8_2()
{
	double a = 0;
	MPI_Init(NULL, NULL);
	int Size; MPI_Comm_size(MPI_COMM_WORLD, &Size);
	int Rank; MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	if (Rank == 0)
	{
		a = 0.3; //"Широковещательная" рассылка сообщения
	}
	MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	if (Rank == 0)
	{
		double* A = new double[Size];
		MPI_Gather(&a, 1, MPI_DOUBLE, A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int k = 0; k < Size; k++)
			cout << A[k] << endl;
		delete[] A;
	}
	else
	{
		double b = a + Rank;
		MPI_Gather(&b, 1, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();
}

void UsrFunc(_In_count_(*len) void* In, _Inout_ void* InOut, _In_ int* len, _In_ MPI_Datatype* Ty)
{
	int* InV = (int*)In;
	int* InOutV = (int*)InOut;
	for (int k = 0; k < (*len); k++)
		InOutV[k] = (InV[k] + InOutV[k]) % 8;
}

void task_8_4()
{
	int N = 8;
	int* A = new int[N];
	MPI_Init(NULL, NULL);
	//int Size=MPI::COMM_WORLD.Get_size();
	int Rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	for (int k=0; k<N; k++)
		A[k]=k+Rank;
	MPI_Op My_Op;
	MPI_Op_create(&UsrFunc, true, &My_Op);
	if (Rank == 0)
	{
		int* B = new int[N];
		MPI_Reduce(A, B, N, MPI_INT, My_Op, 0, MPI_COMM_WORLD);
		for (int k = 0; k < N; k++)
			cout << B[k] << " ";
		cout << endl;
		delete[] B;
	}
	else
	{
		MPI_Reduce(A, NULL, N, MPI_INT, My_Op, 0, MPI_COMM_WORLD);
	}
	MPI_Op_free(&My_Op);
	MPI_Finalize();
	delete[] A;
}

void task_8_6()
{
	MPI_Init(NULL, NULL);
	int Size;
	MPI_Comm_size(MPI_COMM_WORLD, &Size);
	int* Index = new int[Size];
	for (int k = 0; k < Size; k++)
		Index[k] = Size - 1 + k;
	int* Edges = new int[2 * Size - 2];
	for (int k = 0; k < Size - 1; k++)
	{
		Edges[k] = k + 1;
		Edges[Size - 1 + k] = 0;
	}
	MPI_Comm Gr_Comm;
	MPI_Graph_create(MPI_COMM_WORLD, Size, Index, Edges, true, &Gr_Comm);
	int Rank;
	MPI_Comm_rank(Gr_Comm, &Rank); 
	int* Neighbors = new int[Size - 1];
	int Num;
	MPI_Graph_neighbors_count(Gr_Comm, Rank, &Num);
	MPI_Graph_neighbors(Gr_Comm, Rank, Num, Neighbors);
	for (int k = 0; k < Num; k++)
	{
		int Rank1;
		MPI_Status St;
		MPI_Sendrecv(&Rank, 1, MPI_INT, Neighbors[k], 
			5, &Rank1, 1, MPI_INT, Neighbors[k], 5, Gr_Comm, &St);
		cout << "Process " << Rank << " co-operates with process " << Rank1 << endl;
	}
	MPI_Comm_free(&Gr_Comm);
	delete[] Neighbors;
	delete[] Edges;
	delete[] Index;
	MPI_Finalize();
}

