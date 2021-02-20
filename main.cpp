#include <iostream>
#include <math.h>
using namespace std;
#define N  500
#define E 0.05
#define T 0.001


void Printer (double* matrix_A, int n ){
    for (int i = 0; i < n; ++i){
        for (int j = 0; j < n; ++j)
            cout << *(matrix_A + i  * n + j)<< " ";
        cout  << endl;
    }
}


double* Multiplyer_NNxN (const double *A, const double *B, double *C){
    for (int i = 0; i < N; ++i) {
        *(C + i) = 0;
        for (int j = 0; j < N; ++j)
            *(C + i) += *(A + j + N * i) * *(B + j);
    }
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


int main() {
    double matrix_A [N][N];
    double B [N];
    double X_o [N];
    double ptr [N];

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            if (i == j)
                matrix_A[i][j] = 2;
            else
                matrix_A[i][j] = 1;
    }
    for (int i = 0; i < N; ++i) {
        B[i] = N + 1;
        X_o[i] = 0;
    }
   while (chek(Multiplyer_NNxN(&matrix_A[0][0], &X_o[0], &ptr[0]), &B[0]))
       iter (&X_o[0], &ptr[0], &B[0]);

    for (double i : X_o) {
        cout << i <<" ";
    }

    return 0;
}
