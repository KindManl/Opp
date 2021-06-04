
#include <iostream>
#include <cmath>
#include <time.h>
using namespace std;

void MultyMatrix (const int * A, const int * B, int * C, int a, int b, int c){
    //for (int k = 0; k < a; ++k)
    for (int i = 0; i < c; ++i)
        for (int j = 0; j < a; ++j) {
            C[i] += A[j] * B[j * (c - 1) + i];
            //cout << "\n " << A[j] << " " << B[j * (c - 1) + i] << endl;
        }

}


void fill_matrix(int * matrix, int A, int B){
    srand(time(NULL));
    for (int i = 0; i < A * B; i++)
        matrix[i] = 1 + rand() % 5;
}


void PrintMatrix (int * matrix, int A, int B){
    for (int i = 0; i < B; ++i){
        cout << endl;
        for (int j = 0; j < A; ++j)
            cout << matrix[A * i + j] << " ";
    }
    cout << endl;
}


int main() {
    int *matrix_A = NULL;
    int *matrix_B = NULL;
    int *matrix_C = NULL;
    int A = 5;
    int B = 3;
    int C = 4;
    int D = 5;
    matrix_A = new int[A * B];//4 5
    matrix_B = new int[D * C];//5 6
    matrix_C = new int[D * B];
    fill(matrix_C, matrix_C + A * D, 0);
    fill_matrix(matrix_A, A, B);
    fill_matrix(matrix_B, D, C);
    PrintMatrix(matrix_A, A, B);
    PrintMatrix(matrix_B, D, C);
    PrintMatrix(matrix_C, D, B);
    MultyMatrix(matrix_A, matrix_B, matrix_C, A, B, C);
    int b = 0;
    for (int i = 0; i < B; ++i){
        //cout << "\n a = " << matrix_A[i] << " b = " << matrix_B[C * i] << endl;
        b += matrix_A[i] * matrix_B[C * i + 1];
    }
    cout << "\nc[0 0] = " << b << endl;
    PrintMatrix(matrix_C, D, B);
    return 0;
}