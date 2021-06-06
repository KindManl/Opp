#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <mpi/mpi.h>
#include <math.h>

#define W 10 // weight parameter
#define NUM_LISTS 3 // number of lists
#define NUM_TASKS 100   // number of tasks per process
#define NUM_TASKS_TO_SHARE 10

//for threads
pthread_t threads[2];  // thread-descriptors: recv & exec
pthread_mutex_t mutex; // mutex


//global vars
int *tasks; // tasks array
int ProcNum;   // number of processes
int ProcRank;   // rank of current process
int RemainingTasks;  // number of remaining tasks
int ExecutedTasks;   // number of tasks executed
double GlobalRes = 0;   // global result for thread
double GlobalResSumm = 0; // global global result
int STOP_RECV_MARKER = -1; // stop key


void InitTasks(int *tasksArray, int numTasks, int iterCounter) {//non random weight of task now
    for (int i = 0; i < numTasks; i++) {
        tasksArray[i] = abs(50 - i % 100) * abs(ProcRank - (iterCounter % ProcNum)) * W;
    }
}


void ExecuteTasks(const int *tasksArray) {
    for (int i = 0; i < RemainingTasks; i++) {

        // block & unlock task for the proc
        pthread_mutex_lock(&mutex);
        int currTaskWeight = tasksArray[i];
        pthread_mutex_unlock(&mutex);

        for (int j = 0; j < currTaskWeight; j++) {
            GlobalRes += sin(j);
            usleep(tasksArray[i]);  // sleeping task-weight microseconds make weight of task
        }

        ExecutedTasks++;
    }
    RemainingTasks = 0;
}


void* ExecuteRoutine(void* args){
    tasks = (int*)malloc(sizeof(int) * NUM_TASKS);
    MPI_Status status;

    for (int i = 0; i < NUM_LISTS; i++) {

        InitTasks(tasks, NUM_TASKS, i);

        RemainingTasks = NUM_TASKS;
        ExecutedTasks = 0;
        int ExtraTasks = 0;


        ExecuteTasks(tasks);//own tasks

        for (int j = 0; j < ProcNum; j++) {//balance
            if (j != ProcRank) {//not this proc
                MPI_Send(&ProcRank, 1, MPI_INT, j, 0, MPI_COMM_WORLD);

                // extra tasks
                MPI_Recv(&ExtraTasks, 1, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

                if (ExtraTasks) {

                    MPI_Recv(tasks, ExtraTasks, MPI_INT, j, 1, MPI_COMM_WORLD, &status);

                    // executing extra tasks
                    RemainingTasks = ExtraTasks;
                    ExecuteTasks(tasks);
                }
            }
        }

        MPI_Allreduce(&GlobalRes, &GlobalResSumm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Send(&STOP_RECV_MARKER, 1, MPI_INT, ProcRank, 0, MPI_COMM_WORLD);
    free(tasks);
    return NULL;
}

void* ReceiveRoutine(void* args){
    int ShareTasks, requestRank;
    MPI_Status status;

    while (true) {
        //rec request. check request. mute. check tasks. share. unmute. done
        MPI_Recv(&requestRank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
        if (requestRank == STOP_RECV_MARKER) {
            pthread_exit(NULL);
        }


        pthread_mutex_lock(&mutex);
        if (RemainingTasks > NUM_TASKS_TO_SHARE) {
            ShareTasks = RemainingTasks / 2;
            RemainingTasks -= ShareTasks;

            MPI_Send(&ShareTasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
            MPI_Send(&tasks[ShareTasks - 1], ShareTasks, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        } else {
            ShareTasks = 0;
            MPI_Send(&ShareTasks, 1, MPI_INT, status.MPI_SOURCE, 1, MPI_COMM_WORLD);
        }
        pthread_mutex_unlock(&mutex);
    }
}


int CreateThreads() {
    pthread_attr_t attrs;
    //create threads and their attributes
    if (pthread_attr_init(&attrs)){
        perror("COULD NOT INITIALIZE ATTRIBUTES");
        return 1;
    }

    if (pthread_attr_setdetachstate(&attrs, PTHREAD_CREATE_JOINABLE)){
        perror("COULD NOT SET DETACH STATE");
        return 1;
    }

    if (pthread_create(&threads[0], &attrs, ReceiveRoutine, NULL)){//receive thread
        perror("COULD NOT CREATE RECEIVE-THREAD");
        return 1;
    }

    if (pthread_create(&threads[1], &attrs, ExecuteRoutine, NULL)){//execute thread
        perror("COULD NOT CREATE EXECUTE-THREAD");
        return 1;
    }

    pthread_attr_destroy(&attrs);

    if (pthread_join(threads[0], NULL)) {
        perror("COULD NOT JOIN RECEIVE-THREAD");
        return 1;
    }

    if (pthread_join(threads[1], NULL)){
        perror("COULD NOT JOIN EXECUTE-THREAD");
        return 1;
    }
    return 0;
}


void TasksDistribution() {
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    //mutex
    pthread_mutex_init(&mutex, NULL);
    CreateThreads();
    pthread_mutex_destroy(&mutex);
}


int main(int argc, char **argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Cant make threads\n");
        MPI_Finalize();
        return 1;
    }

    TasksDistribution();

    MPI_Finalize();
}