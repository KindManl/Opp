#include <iostream>
#include <cmath>
#include <mpi/mpi.h>
using namespace std;
#define N  5110
#define E 0.0005
#define T 0.0001


double* Multiplyer_NNxN (const double *A, const double *B, double *C, int NLines){
    for (int i = 0; i < NLines; ++i) {
        *(C + i) = 0.0f;
        for (int j = 0; j < N; ++j)
            *(C + i) += *(A + j + N * i) * *(B + j);
    }

    return C;
}


bool check (const double* A, const double* B, int NLines, int FirstLine){
    double temp[NLines];
    double B_glob = 0;//all B_m;
    double B_m = 0;
    double A_glob = 0;
    double A_M = 0;
    for (int i = 0; i < NLines; ++i)
        temp[i] =  - *(B + i + FirstLine) + *(A + i);
    for (int i = 0; i <  NLines; ++i){
        B_m += B[i + FirstLine] * B[i + FirstLine];
        A_M += temp[i] * temp[i];
    }
    MPI_Allreduce(&B_m,&B_glob,1,MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    B_m = sqrt (B_glob);

    MPI_Allreduce(&A_M, &A_glob,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    A_M = sqrt(A_glob);
    return (A_M / B_m) >= E;
}


void  iter (double *X, const double *ptr, const double *B, int NLines, int FirstLine){
    for (int i = 0; i < NLines; ++i)
        *(X + FirstLine + i) = *(X + FirstLine + i) - T * (*(ptr + i) - *(B + FirstLine + i));


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

    for (int i = 0; i < N; ++i)
        X_o[i] = 0;
    for (int i = FirstLine; i < FirstLine + NLines; ++i) {
        B[i] = N + 1;
        ptr[i - FirstLine] = 0;
    }


    while(check(Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0], NLines), &B[0], NLines, FirstLine)){
        iter(&X_o[0], &ptr[0], &B[0], NLines, FirstLine);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allgatherv(X_o, CountEl[ProcRank], MPI_DOUBLE, X_o, CountEl, shift, MPI_DOUBLE, MPI_COMM_WORLD);
    }

    if (LastProc)
        for (double i : X_o)
            cout << i <<" ";

    MPI_Finalize();
    cout<<endl;
    return 0;
}
