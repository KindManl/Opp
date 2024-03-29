#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
#define N  500 //num of lines in all processes
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
    double SUMM = 0;
    double ptr;
    for (int i = 0; i < Number_Of_Lines; ++i) {
        ptr = A[i] / B[FirstLine + i];
        SUMM += ptr * ptr - 2 * ptr + 1;
    }
    MPI_Allreduce(&SUMM, &SUMM,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return SUMM >= E * E;
}

void  iteration (double *X, const double *ptr, const double *B, int Number_Of_Lines, int FirstLine){//find next X value
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
    fill(CountEl, CountEl + ProcNum - 1, div_up(N, ProcNum));
    CountEl[ProcNum - 1] = N - div_up(N, ProcNum) * (ProcNum - 1);

    //number of first line in matrix from big matrix
    int FirstLine = Number_Of_Lines * ProcRank;
    if (LastProc)
        FirstLine = N - Number_Of_Lines;

    int shift [ProcNum];//displs for allgather
    for (int i = 0; i < ProcNum; ++i)
        shift[i] = div_up(N, ProcNum) * i;

    double *matrix_A = new double [N * Number_Of_Lines];//matrix A
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
    fill(X_o, X_o + N, 0);
    fill(B + FirstLine, B + FirstLine + Number_Of_Lines, N + 1);
    fill(ptr, ptr + Number_Of_Lines, 0);

    //calculations block
    do{
        Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0], Number_Of_Lines);
        iteration(&X_o[0], &ptr[0], &B[0], Number_Of_Lines, FirstLine);
        MPI_Allgatherv(X_o, CountEl[ProcRank], MPI_DOUBLE, X_o, CountEl, shift, MPI_DOUBLE, MPI_COMM_WORLD);
    } while (check(&ptr[0], &B[0], Number_Of_Lines, FirstLine));

    if (LastProc)
        for (int i = 0; i < N; ++i)
            cout <<X_o[i]<<" ";

    MPI_Finalize();
    return 0;
}
