#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
#define N  5110 //num of lines in all processes
#define E 0.0005
#define T 0.0001


double* Multiplyer_NNxN (const double *A, const double *B, double *C, int Number_Of_Lines){
    for (int i = 0; i < Number_Of_Lines; ++i) {
        C[i] = 0.0f;
        for (int j = 0; j < N; ++j)
            C[i] += A[j + N * i] * B[j];
    }

    return C;
}


bool check (const double* A, const double* B, int Number_Of_Lines, int FirstLine){
    double temp[Number_Of_Lines];
    double B_m = 0;
    double A_M = 0;
    for (int i = 0; i < Number_Of_Lines; ++i)
        temp[i] =  - B[i + FirstLine] + A[i];

    for (int i = 0; i < Number_Of_Lines; ++i){
        B_m += B[i + FirstLine] * B[i + FirstLine];
        A_M += temp[i] * temp[i];
    }

    B_m /= Number_Of_Lines;
    B_m = sqrt (B_m);

    A_M /= Number_Of_Lines;
    A_M = sqrt(A_M);
    return (A_M / B_m) >= E;
}


void  iter (double *X, const double *ptr, const double *B, int Number_Of_Lines, int FirstLine){
    for (int i = 0; i < Number_Of_Lines; ++i)
        X[FirstLine + i] = X[FirstLine + i] - T * (ptr[i] - B[i + FirstLine]);

}


int div_up (int a, int b){
    return (a + b - 1) / b;
}


int main(int argc, char *argv[]) {
    //MPI init block
    MPI_Init(&argc, &argv);
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD,&ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD,&ProcRank);

    //last proc block && Num lines in proc
    bool LastProc = false;
    if (ProcRank == ProcNum - 1)
        LastProc = true;
    int Number_Of_Lines = div_up(N, ProcNum);//Num of lines in each proc
    if (LastProc) {//number of lines in proc, last proc less lines then others
        Number_Of_Lines = N - div_up(N, ProcNum) * (ProcNum - 1);
    }

    //buff of number of elements for allgatherv func
    int CountEl [ProcNum];
    int id = div_up(N, ProcNum);
    for (int i = 0; i < ProcNum - 1; ++i)
        CountEl[i]  = id;
    CountEl[ProcNum - 1] = N - div_up(N, ProcNum) * (ProcNum - 1);

    //number of first line in matrix from big matrix
    int FirstLine = Number_Of_Lines * ProcRank;
    if (LastProc)
        FirstLine = N - Number_Of_Lines;


    double *matrix_A = nullptr;
    matrix_A = new double [N * Number_Of_Lines];

    int shift [ProcNum];//displs for allgather
    for (int i = 0; i < ProcNum; ++i)
        shift[i] = 0 + div_up(N, ProcNum) * i;

    double *B = new double [N];//B vector
    double *X_o = new double [N];//result vector in proc & task
    double *ptr = new double [Number_Of_Lines];//ptr for multy and chek func

    //buffers init block
    for (int i = 0; i < Number_Of_Lines; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j - Number_Of_Lines * ProcRank)
                matrix_A[i * N + j] = 2.0f;
            else
                matrix_A[i * N + j] = 1.0f;
        }
    }

    for (int i = 0; i < N; ++i) {
        B[i] = N + 1;
        X_o[i] = 0;
    }
    while  (check(Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0], Number_Of_Lines), &B[0], Number_Of_Lines, FirstLine)) {
        iter(&X_o[0], &ptr[0], &B[0], Number_Of_Lines, FirstLine);
        MPI_Allgatherv(X_o, CountEl[ProcRank], MPI_DOUBLE, X_o, CountEl, shift, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    if (LastProc)
        for (int i = 0; i < N; ++i)
            cout <<X_o[i]<<" ";

    MPI_Finalize();
    cout<<endl;
    return 0;
}
