#include <iostream>
#include <cmath>
#include <omp.h>
using namespace std;
#define N  1000
#define E 0.00001
#define T 0.000001


int div_up (int a, int b){
    return (a + b - 1) / b;
}


int main(){
    double *matrix_A = new double [N * N];//matrix
    double *B = new double [N];           //B vector
    double *X_o = new double [N];         //result vector in proc & task
    double *ptr = new double [N];         //ptr for multy and chek func
    double *tmp = new double [N];         //tmp for norma func
    double Norma_A;
    double Norma_B;
    double Epsilon = E * E;

    #pragma omp parallel for
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

    #pragma omp parallel
    do{
    #pragma omp for
        for (int i = 0; i < N; ++i) {
            ptr[i] = 0.0f;
            for (int j = 0; j < N; ++j)
                ptr[i] += matrix_A[j + N * i] * X_o[j];
        }

        #pragma omp for
        for (int i = 0; i < N; ++i)
            X_o[i] = X_o[i] - T * (ptr[i] - B[i]);

        #pragma omp for
        for (int i = 0; i < N; ++i) {
            ptr[i] = (ptr[i] - B[i]) * (ptr[i] - B[i]);
            tmp[i] = B[i] * B[i];
        }

        Norma_A = 0;
        Norma_B = 0;

        #pragma omp for
        for (int i = 0; i < N; ++i){
            Norma_A += ptr[i] * ptr[i];
            Norma_B += B[i] * B[i];
        }

    } while (Norma_A / Norma_B >= Epsilon);


    for (int i = 0; i < N; ++i)
        cout << X_o[i] <<" ";
    cout<<endl;
    return 0;
}