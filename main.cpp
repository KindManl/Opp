#include <iostream>
#include <cmath>
#include <mpich/mpi.h>
using namespace std;
#define N  500
#define E 0.005
#define T 0.00001



double* Multiplyer_NNxN (const double *A, const double *B, double *C, int ProcRank, int NLines, const int* shift, const int *Count){
    for (int i = 0; i < NLines; ++i) {
        *(C + i) = 0.0f;
        for (int j = 0; j < N; ++j)
            *(C + i) += *(A + j + N * i) * *(B + j);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(C, Count[ProcRank], MPI_DOUBLE, C, Count, (shift), MPI_DOUBLE, MPI_COMM_WORLD);

    return C;
}


bool chek (const double* A, const double* B){
    double temp[N];
    double B_m = 0;
    double A_M = 0;
    for (int i = 0; i < N; ++i)
        temp[i] =  - *(B + i) + *(A + i);

    for (int i = 0; i < N; ++i){
        B_m += B[i] * B[i];
        A_M += temp[i] * temp[i];
    }

    B_m /= N;
    B_m = sqrt (B_m);

    A_M /= N;
    A_M = sqrt(A_M);
    return (A_M / B_m) >= E;
}

void  iter (double *X, const double *ptr, const double *B){
    for (int i = 0; i < N; ++i)
        *(X + i) = *(X + i) - T * (*(ptr + i) - *(B + i));

}


int div_up (int a, int b){
    return (a + b - 1) / b;
}


int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank);

    bool LastProc;
    if (ProcRank == ProcNum)
        LastProc = true;
    int NLines = div_up(N, ProcNum);
    if (LastProc) {//number of lines in proc, last proc less lines then others
        NLines = N - div_up(N, ProcNum) * (ProcNum);
    }

    int CountEl [ProcNum];
    int id = div_up(N, ProcNum);
    for (int i = 0; i < ProcNum - 1; ++i)
        CountEl[i]  = id;
    CountEl[ProcNum] = N - div_up(N, ProcNum) * (ProcNum);


    int FirstLine = NLines * ProcRank;//start line in big matrix
    if (LastProc)
        FirstLine = N - NLines;

    double *matrix_A = nullptr;
    matrix_A = new double [N * NLines];
    int shift = div_up(N, ProcNum);
    double B [N];
    double X_o [N];
    double ptr [N];
    for (int i = 0; i < NLines; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j + NLines * ProcRank)
                matrix_A[i * N + j] = 2.0f;
            else
                matrix_A[i * N + j] = 1.0f;
        }
    }

    for (int i = 0; i < N; ++i) {
        B[i] = N + 1;
        X_o[i] = 0;
    }


    while (chek(Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0], ProcRank, NLines, &shift, &CountEl[0]), &B[0]))
        iter (&X_o[0], &ptr[0], &B[0]);

    for (double i : X_o) {
        cout << i <<" ";
    }
    MPI_Finalize();
    return 0;
}
