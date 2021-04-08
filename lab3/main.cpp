#include <iostream>
#include <mpi/mpi.h>
#include <cmath>
#include <time.h>
using namespace std;


bool check_size(int A, int B, int C, int D){
    if ((A < 0) || (B < 0) || (C < 0) || (D < 0)){
        cout << "\nMatrix size error\n";
        return false;
    }
    if (B != C){
        cout << "\nB should be equal C\n";
        return false;
    }
    return true;
}


bool AllocData (int* A, int* B, int* C, int* D){
    cout << "enter size of matrix 1" << endl;
    cin >> *A >> *B;
    cout << "enter size of matrix 2" << endl;
    cin >> *C >> *D;
    return check_size(*A, *B, *C, *D);
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
    cout << "h";
    for (int i = 0; i < a * b; ++i)
        B[i] = A[i];
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

    MPI_Comm_size (MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank (MPI_COMM_WORLD, &ProcRank);

    MPI_Dims_create(ProcNum, 2, dims);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &BaseComm);
    MPI_Cart_coords(BaseComm, ProcRank, 2, coords);
    MPI_Comm_split(BaseComm, coords[1], coords[0], &rowComm);
    MPI_Comm_split(BaseComm, coords[0], coords[1], &colComm); //make smthnk with this shit

    int *matrix_A = NULL;
    int *matrix_B = NULL;
    int *matrix_C = NULL;
    int *segment_A = NULL;
    int *segment_B = NULL;
    int *segment_C = NULL;
    int A, B, C, D, Rows, Columns, tmp [3];

    if (ProcRank == 0) {
        if (!AllocData(&A, &B, &C, &D)){
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 0;
        }

        matrix_A = new int[A * B];
        matrix_B = new int[C * D];
        matrix_C = new int[A * D];

        fill(matrix_C, matrix_C + A * D, 0);
        fill_matrix(matrix_A, A, B);
        fill_matrix(matrix_B, C, D);

        tmp[0] = A;
        tmp[1] = D;
        tmp[2] = B;
    }

    MPI_Bcast(tmp, 3, MPI_INT, 0, MPI_COMM_WORLD);//sharing A B C values with all procs
    A = tmp[0];
    D = tmp[1];
    B = tmp[2];

    Rows = A / dims[0];//Rows & columns of segments in procs
    Columns = D / dims[1];

    cout << "\n rows = " << Rows << " col = " << Columns << " 0 = " << dims[0] << " 1 = " << dims[1] << " r = " << ProcRank << endl;
    segment_A = new int[Rows * B];
    segment_B = new int[B * Columns];
    segment_C = new int[Rows * Columns];

    //scate matrix A with rows
    MPI_Scatter(matrix_A, A * B / dims[0], MPI_INT, segment_A, B * Rows, MPI_INT, 0, rowComm);
    MPI_Bcast(segment_A, Rows * B, MPI_INT, 0, colComm);

    //scate matrix B with columns
    MPI_Datatype COLUMN;
    MPI_Type_vector(B, 1, D, MPI_INT, &COLUMN);
    MPI_Type_commit(&COLUMN);
    // MPI_Scatter(matrix_B, 1, COLUMN, segment_B, Columns * B, MPI_INT, 0, colComm);

    MPI_Type_free(&COLUMN);


    if (ProcRank == 1)
        PrintMatrix(segment_A, Rows, B);
    if (ProcRank == 0)
        PrintMatrix(segment_A, Rows, B);

    //MultyMatrix(segment_A, segment_B, segment_C, Rows, B, Columns);
    if (ProcRank == 0) {
        
        cout << "\npro 0" << endl;
        PrintMatrix(matrix_A, A, B);
        //PrintMatrix(matrix_B, C, D);
        //MultyMatrix(matrix_A, matrix_B, matrix_C, A, B, D);
        //PrintMatrix(segment_C, Rows, Columns);
        PrintMatrix(segment_A, Rows, B);
    }
    if (ProcRank == 1) {

        cout << "\npro 1" << endl;
        //PrintMatrix(matrix_A, A, B);
        //PrintMatrix(matrix_B, C, D);
        //MultyMatrix(matrix_A, matrix_B, matrix_C, A, B, D);
        //PrintMatrix(segment_C, Rows, Columns);
        PrintMatrix(segment_A, Rows, B);
    }
    MPI_Finalize();
    return 0;
}
