#include <mpi/mpi.h>
#include <iostream>
#include <cmath>


#define A 0.00001
#define E 0.00000001
#define Nodes 64

double finding_func(double x, double y, double z) {//finding func
    //work
    return x * x + y * y + z * z;
}


double ro(double x, double y, double z) {//calculates ro from task
    //work
    return 6 - A * finding_func(x, y, z);
}


double finding_func_inside(const int x, const int y, const int z, double* cube, const double incremention, const int ProcNum, const int ProcRank){//calculates values of finding func inside the cube
    //must work
    double x_start = ProcRank * Nodes / ProcNum;
    int a, b, c, d, e, f;
    a = ((x - 1) * Nodes * Nodes) + (y * Nodes) + z;
    b = ((x + 1) * Nodes * Nodes) + (y * Nodes) + z;
    c = (x * Nodes * Nodes) + ((y - 1) * Nodes) + z;
    d = (x * Nodes * Nodes) + ((y + 1) * Nodes) + z;
    e = (x * Nodes * Nodes) + (y * Nodes) + z - 1;
    f = (x * Nodes * Nodes) + (y * Nodes) + z + 1;

    double BoardSumm = cube[a] + cube[b] + cube[c] + cube[d] + cube[e] + cube[f];
    double ptr = BoardSumm / incremention + ro((-1.0 + (x_start + x)) * incremention, -1.0 + y * incremention, -1.0 + z * incremention);
    ptr = ptr * 6 / (incremention * incremention + A);
    return 1.0 / ptr;
}

bool check(){
    return false;
}



int main(int argc, char* argv[])
{
    //loc vars for MPI and calculations
    double start, end, check_value, diff_value, max_diff_value;
    int local_stop_flag = 0;
    int global_stop_flag = 0;
    int ProcNum, ProcRank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Request snd_request[2];
    MPI_Request rcv_request[2];


    //set cube
    double* cube = (double*)calloc((Nodes / ProcNum + 2) * Nodes * Nodes, sizeof(double));
    double* cube_old = (double*)calloc((Nodes / ProcNum + 2) * Nodes * Nodes, sizeof(double));
    int x_starting_coord = ProcRank * Nodes / ProcNum;
    double incremention = 2.0 / Nodes;
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


    //set old cube
    for (int i = 0; i < (Nodes / ProcNum + 2) * Nodes * Nodes; ++i)
        cube_old[i] = 0;


    while(true){
    //for (int i = 0; i < 4; ++i){
        max_diff_value = 0;
        //send receive bloc
        if (ProcRank != 0) {
            MPI_Isend(cube + (Nodes * Nodes), Nodes * Nodes, MPI_DOUBLE, ProcRank - 1, 0, MPI_COMM_WORLD, &snd_request[0]);
            MPI_Irecv(cube, Nodes * Nodes, MPI_DOUBLE, ProcRank - 1, 1, MPI_COMM_WORLD, &rcv_request[1]);
        }
        if (ProcRank != ProcNum - 1) {
            MPI_Isend(cube + (ProcNum / Nodes * Nodes * Nodes), Nodes * Nodes, MPI_DOUBLE, ProcRank + 1, 1, MPI_COMM_WORLD, &snd_request[1]);
            MPI_Irecv(cube + (Nodes * Nodes * (ProcNum / Nodes + 1)), Nodes * Nodes, MPI_DOUBLE, ProcRank + 1, 0, MPI_COMM_WORLD, &rcv_request[0]);
        }

        //calculate bloc
        for (int x = (Nodes / ProcNum + 2) / 2; x > 1; x--) {
            for (int y = 1; y < Nodes; y++) {
                for (int z = 1; z < Nodes; z++) {
                    cube[x * Nodes * Nodes + y * Nodes + z] = finding_func_inside(x, y, z, cube, incremention, ProcNum, ProcRank);
                }
            }
        }
        for (int x = (Nodes / ProcNum + 2) / 2; x < Nodes / ProcNum; x++) {
            for (int y = 1; y < Nodes; y++) {
                for (int z = 1; z < Nodes; z++) {
                    cube[x * Nodes * Nodes + y * Nodes + z] = finding_func_inside(x, y, z, cube, incremention, ProcNum, ProcRank);
                }
            }
        }


        //wait bloc
        if (ProcRank != 0) {
            MPI_Wait(&snd_request[0], MPI_STATUS_IGNORE);
            MPI_Wait(&rcv_request[1], MPI_STATUS_IGNORE);
        }
        if (ProcRank != ProcNum - 1) {
            MPI_Wait(&snd_request[1], MPI_STATUS_IGNORE);
            MPI_Wait(&rcv_request[0], MPI_STATUS_IGNORE);
        }



        //calculate bloc 2 (boarders)
        if (ProcRank != 0) {
            for (int i = 1; i < Nodes - 1; i++) {
                for (int j = 1; j < Nodes - 1; j++) {
                    cube[(1 * Nodes * Nodes) + (Nodes * j) + i] = finding_func_inside(1, j, i, cube, incremention, ProcNum, ProcRank);
                }
            }
        }
        if (ProcRank != ProcNum - 1) {
            for (int z = 1; z < Nodes - 1; z++) {
                for (int y = 1; y < Nodes - 1; y++) {
                    cube[((Nodes / ProcNum) * Nodes * Nodes) + (Nodes * y) + z] = finding_func_inside(Nodes / ProcNum, y, z, cube, incremention, ProcNum, ProcRank);
                }
            }
        }


        //check bloc
        for (int i = 0; i < (Nodes / ProcNum + 2) * Nodes * Nodes; ++i){//need borders to check
            if (cube[i] != 0) {
                diff_value = fabs(cube_old[i] - cube[i]);
                if (diff_value > max_diff_value)
                    max_diff_value = diff_value;
            }
        }

        if (max_diff_value < E) {
            local_stop_flag = 0;
        } else {
            local_stop_flag = 1;
        }
        MPI_Allreduce(&local_stop_flag, &global_stop_flag, 1, MPI_INT, MPI_MAX,
                      MPI_COMM_WORLD);
        if (global_stop_flag == 0) {
            break;
        }


        //reset old cube
        for (int i = 0; i < (Nodes / ProcNum + 2) * Nodes * Nodes; ++i)
            cube_old[i] = cube[i];
    }


    free(cube);
    free(cube_old);
    MPI_Finalize();
    return 0;
}
