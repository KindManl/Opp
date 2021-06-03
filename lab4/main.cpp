#include <mpi/mpi.h>
#include <iostream>


#define a 1e5
#define Nodes 64

double finding_func(double x, double y, double z) {//finding func
    return x * x + y * y + z * z;
}


double ro(double x, double y, double z) {//calculates ro from task
    return 6 - a * finding_func(x, y, z);
}


double finding_func_inside(){//calculates values of finding func inside the cube
    double ptr = 4;
    return ptr;
}

bool check(){
    return false;
}


int main(int argc, char* argv[])
{
    //loc vars for MPI and calculations
    double start, end, check_value;
    int ProcNum, ProcRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    //set cube
    double* cube = (double*)calloc((Nodes / ProcNum + 2) * Nodes * Nodes, sizeof(double));
    int x_starting_coord = ProcRank * Nodes / ProcNum;
    double incremention = 2.0 / Nodes;
    double borders_x[2] = { -1.0, 1.0 };
    double borders_y[2] = { -1.0, 1.0 };
    double borders_z[2] = { -1.0, 1.0 };

    //set boundary values
    for (int z = 0; z < 2; z++) {
        for (int x = 1; x < Nodes / ProcNum + 1; x++) {
            for (int y = 0; y < Nodes; y++) {
                cube[(x * Nodes * Nodes) + (Nodes * (z % 2)) + (y * Nodes)] = finding_func(-1.0 + incremention * (x_starting_coord + x), -1.0 + incremention * y, borders_z[z]);
            }
        }
    }
    for (int y = 0; y < 2; y++) {
        for (int x = 1; x < Nodes / ProcNum + 1; x++) {
            for (int z = 0; z < Nodes; z++) {
                cube[(x * Nodes * Nodes) + ((Nodes - 1) * (y % 2)) + z] = finding_func(-1.0 + incremention * (x_starting_coord + x), borders_y[y], -1.0 + incremention * z);
            }
        }
    }
    if (ProcRank == 0) {
        for (int z = 0; z < Nodes; z++) {
            for (int y = 0; y < Nodes; y++) {
                cube[(1 * Nodes * Nodes) + (Nodes * y) + z] = finding_func(-1, -1.0 + incremention * y, -1.0 + incremention * y);
            }
        }
    } else if (ProcRank == ProcNum - 1) {
        for (int z = 0; z < Nodes; z++) {
            for (int y = 0; y < Nodes; y++) {
                cube[((Nodes / ProcNum) * Nodes * Nodes) + (Nodes * y) + z] = finding_func(1, -1.0 + incremention * y, -1.0 + incremention * y);
            }
        }
    }

    while(check()){
        finding_func_inside();
    }


    free(cube);
    MPI_Finalize();
    return 0;
}
