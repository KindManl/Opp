#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
#define N  527
#define E 0.0005
#define T 0.0001


double* Multiplyer_NNxN (const double *A, const double *B, double *C, int NLines, int FirstLine){
    for (int i = 0; i < NLines; ++i) {
        *(C + i) = 0.0f;
        for (int j = 0; j < N; ++j)
            *(C + i) += *(A + j + N * i) * *(B + i + FirstLine);
    }

    return C;
}


bool chek (const double* A, const double* B, int NLines, int FirstLine){
    double temp[NLines];
    double B_m = 0;
    double A_M = 0;
    for (int i = 0; i < NLines; ++i)
        temp[i] =  - *(B + i + FirstLine) + *(A + i);

    for (int i = 0; i <  NLines; ++i){
        B_m += B[i + FirstLine] * B[i + FirstLine];
        A_M += temp[i] * temp[i];
    }

    B_m /= NLines;
    B_m = sqrt (B_m);

    A_M /= NLines;
    A_M = sqrt(A_M);
    return (A_M / B_m) >= E;
}


void  iter (double *X, const double *ptr, const double *B, int NLines, int FirstLine){
    for (int i = 0; i < NLines; ++i)
        *(X + FirstLine + i) = *(X + FirstLine + i) - T * (*(ptr + i) - *(B + i));


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
    int NLines = div_up(N, ProcNum);
    if (LastProc) {//number of lines in proc, last proc less lines then others
        NLines = N - div_up(N, ProcNum) * (ProcNum - 1);
    }

    //buff of number of elements for allgatherv func
    int CountEl [ProcNum];
    int id = div_up(N, ProcNum);
    for (int i = 0; i < ProcNum - 1; ++i)
        CountEl[i]  = id;
    CountEl[ProcNum - 1] = N - div_up(N, ProcNum) * (ProcNum - 1);

    //number of first line in matrix from big matrix
    int FirstLine = NLines * ProcRank;
    if (LastProc)
        FirstLine = N - NLines;


    double *matrix_A = nullptr;
    matrix_A = new double [N * NLines];

    int shift [ProcNum];//displs for allgather
    for (int i = 0; i < ProcNum; ++i)
        shift[i] = 0 + div_up(N, ProcNum) * i;


    double B [N];//B vector
    double X_o [N];//result vector in proc
    double X [N];//result vector in task
    double ptr [NLines];//ptr for multy and chek func

    //buffers init block
    for (int i = 0; i < NLines; ++i) {
        for (int j = 0; j < N; ++j) {
            if (i == j - NLines * ProcRank)
                matrix_A[i * N + j] = 2.0f;
            else
                matrix_A[i * N + j] = 1.0f;
        }
    }

    for (int i = 0; i < N; ++i) {
        B[i] = N + 1;
        X_o[i] = 0;
    }

    while (chek(Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0], NLines, FirstLine), &B[0], NLines, FirstLine))
        iter(&X_o[0], &ptr[0], &B[0], NLines, FirstLine);


    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgatherv(X_o, CountEl[ProcRank], MPI_DOUBLE, X_o, CountEl, shift, MPI_DOUBLE, MPI_COMM_WORLD);

    if (LastProc)
        for (double i : X_o)
            cout << i <<" ";
        
    MPI_Finalize();
    cout<<endl;
    return 0;
}
