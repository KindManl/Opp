#include <iostream>
#include <mpi/mpi.h>
#include <cmath>
#include <time.h>
using namespace std;


bool check_size(int A, int B, int C, int D, int ProcNum){
    if ((A < 0) || (B < 0) || (C < 0) || (D < 0)){
        cout << "\nMatrix size error\n";
        return false;
    }
    if ((A % ProcNum != 0) || (D % ProcNum)){
        cout << "\n A and D should devide on procnum";
        return false;
    }
    if (B != C){
        cout << "\nB should be equal C\n";
        return false;
    }
    return true;
}


bool AllocData (int* A, int* B, int* C, int* D, int ProcNum){
    cout << "enter size of matrix 1" << endl;
    cin >> *A >> *B;
    cout << "enter size of matrix 2" << endl;
    cin >> *C >> *D;
    return check_size(*A, *B, *C, *D, ProcNum);
}


void fill_matrix(int * matrix, int A, int B){
    srand(time(NULL));
    for (int i = 0; i < A * B; i++)
        matrix[i] = 1 + rand() % 3;
}


void PrintMatrix (int * matrix, int A, int B){
    for (int i = 0; i < A; ++i){
        cout << endl;
        for (int j = 0; j < B; ++j)
            cout << matrix[B * i + j] << " ";
    }
    cout << endl;
}


void MultyMatrix (const int * A, const int * B, int * C, int a, int b, int c){
    for (int k = 0; k < a; ++k)
        for (int i = 0; i < c; ++i)
            for (int j = 0; j < b; ++j)
                C[k * c + i] += A[k * b + j] * B[j * c + i];
}


void CopyMatrix (const int * A, int * B, int a, int b){
    for (int i = 0; i < a * b; ++i)
        B[i] = A[i];
}


void Transposition (int * A, int a, int b){//a 5 b 6
    int * tmp = new int [a * b];
    for (int i = 0; i < b; ++i)
        for (int j = 0; j < a; ++j)
            tmp[i * a + j] = A[i + j * b];
    CopyMatrix(tmp, A, a, b);
    delete[] tmp;
}


int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm BaseComm;
    MPI_Comm rowComm;
    MPI_Comm colComm;

    int dims[2] = {0, 0};
    int periods[2] = {0, 0};
    int coords[2];
    int reorder = 0;
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    MPI_Dims_create(ProcNum, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &BaseComm);
    MPI_Cart_coords(BaseComm, ProcRank, 2, coords);
    MPI_Comm_split(BaseComm, coords[1], coords[0], &rowComm);
    MPI_Comm_split(BaseComm, coords[0], coords[1], &colComm);

    int *matrix_A = NULL;
    int *matrix_B = NULL;
    int *matrix_C = NULL;
    int *segment_A = NULL;
    int *segment_B = NULL;
    int *segment_C = NULL;
    int A, B, C, D, Rows, Columns, tmp[3];

    if (ProcRank == 0) {
        if (!AllocData(&A, &B, &C, &D, ProcNum)) {
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        matrix_A = new int[A * B];
        matrix_B = new int[C * D];
        matrix_C = new int[A * D];

        fill(matrix_C, matrix_C + A * D, 0);
        fill_matrix(matrix_A, A, B);
        fill_matrix(matrix_B, C, D);

        cout << "\nMatrix A:";
        PrintMatrix(matrix_A, A, B);
        cout << "\nMatrix B:";
        PrintMatrix(matrix_B, C, D);
        Transposition(matrix_B, C, D);

        tmp[0] = A;
        tmp[1] = D;
        tmp[2] = B;
    }

    MPI_Bcast(tmp, 3, MPI_INT, 0, MPI_COMM_WORLD);//sharing A B C values with all procs
    A = tmp[0];
    D = tmp[1];
    B = tmp[2];

    Rows = A / dims[1];//Rows & columns of segments in procs
    Columns = D / dims[0];

    segment_A = new int[Rows * B];
    segment_B = new int[B * Columns];
    segment_C = new int[Rows * Columns];
    fill(segment_C, segment_C + Rows * Columns, 0);

    //scate matrix A with rows
    if (coords[0] == 0)
        MPI_Scatter(matrix_A, Rows * B, MPI_INT, segment_A, B * Rows, MPI_INT, 0, colComm);
    MPI_Bcast(segment_A, Rows * B, MPI_INT, 0, rowComm);

    //scate matrix B with columns
    if (coords[1] == 0) {
        MPI_Scatter(matrix_B, B * Columns, MPI_INT, segment_B, B * Columns, MPI_INT, 0, rowComm);
        Transposition(segment_B, Columns, B);
    }
    MPI_Bcast(segment_B, B * Columns, MPI_INT, 0, colComm);

    MultyMatrix(segment_A, segment_B, segment_C, Rows, B, Columns);

    MPI_Datatype Receiver;
    MPI_Datatype Receiver_elem;
    MPI_Type_vector(Rows, Columns, D, MPI_INT, &Receiver_elem);
    MPI_Type_create_resized(Receiver_elem, 0, Columns * sizeof(int), &Receiver);
    MPI_Type_commit(&Receiver);

    MPI_Gather(segment_C, Rows * Columns, MPI_INT, matrix_C, 1, Receiver, 0, BaseComm);

    MPI_Type_free(&Receiver);
    MPI_Type_free(&Receiver_elem);

    if (ProcRank == 0) {
        cout << "\nMatrix C:";
        PrintMatrix(matrix_C, A, D);
    }
    delete[] matrix_A;
    delete[] matrix_B;
    delete[] matrix_C;
    delete[] segment_A;
    delete[] segment_B;
    delete[] segment_C;
    MPI_Comm_free(&rowComm);
    MPI_Comm_free(&colComm);
    MPI_Comm_free(&BaseComm);
    MPI_Finalize();
    return 0;
}
