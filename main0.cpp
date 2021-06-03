#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;
#define N  500
#define E 0.00001
#define T 0.000001


double* Multiplyer_NNxN (const double *A, const double *B, double *C){
    for (int i = 0; i < N; ++i) {
        C[i] = 0.0f;
        for (int j = 0; j < N; ++j)
            C[i] += A[j + N * i] * B[j];
    }
    return C;
}


bool check (const double* A, const double* B){
    double SUMM = 0;
    double ptr;
    for (int i = 0; i < N; ++i) {
        ptr = A[i] / B[i];
        SUMM += ptr * ptr - 2 * ptr + 1;
    }
    return SUMM >= E * E;
}


void  iteration (double *X, const double *ptr, const double *B){
    for (int i = 0; i < N; ++i)
        X[i] = X[i] - T * (ptr[i] - B[i]);
}


int div_up (int a, int b){
    return (a + b - 1) / b;
}


int main(){
    double *matrix_A = new double [N * N];//matrix
    double *B = new double [N];           //B vector
    double *X_o = new double [N];         //result vector in proc & task
    double *ptr = new double [N];         //ptr for multy and chek func

//#pragma omp parallel for
    for (int i = 0; i < N; ++i) {//init block
        for (int j = 0; j < N; ++j) {
            if (i == j)
                matrix_A[i * N + j] = 2.0f;
            else
                matrix_A[i * N + j] = 1.0f;
        }
        B[i] = N + 1;
        X_o[i] = 0;
    }
    int k = 0;
    while (check(Multiplyer_NNxN(&matrix_A[0], &X_o[0], &ptr[0]), &B[0])) {
        iteration(&X_o[0], &ptr[0], &B[0]);
        k++;
    }

    for (int i = 0; i < N; ++i)
        cout << X_o[i] <<" ";
    cout<<endl;
    return 0;
}