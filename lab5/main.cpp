#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>
#include <mpi/mpi.h>
#include <math.h>
#include <stdbool.h>

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
double GlobalRes = 0;   // global result


void InitTasks(int *tasksArray, int numTasks, int iterCounter) {//non random weight of task now
    for (int i = 0; i < numTasks; i++) {
        tasksArray[i] = abs(50 - i % 100) * abs(rank - (iterCounter % numProcesses)) * W;
    }
}


void ExecuteTasks(const int *tasksArray) {
    for (int i = 0; i < RemainingTasks; i++) {

        // block & unlock task for the proc
        pthread_mutex_lock(&mutex);
        int currTaskWeight = tasksArray[i];
        pthread_mutex_unlock(&mutex);

        for (int j = 0; j < currTaskWeight; j++) {
            globalRes += sin(j);
            usleep(tasksArray[i]);  // sleeping task-weight microseconds make weight of task
        }

        ExecutedTasks++;
    }
    numRemainingTasks = 0;
}


void TasksDistribution() {
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //mutex
    pthread_mutex_init(&mutex, NULL);
    createThreads();
    pthread_mutex_destroy(&mutex);
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

    if (pthread_create(&threads[0], &attrs, receiveStartRoutine, NULL)){//receive thread
        perror("COULD NOT CREATE RECEIVE-THREAD");
        return 1;
    }

    if (pthread_create(&threads[1], &attrs, executeStartRoutine, NULL)){//execute thread
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

void* ExecuteRoutine(void* args){
    //work
}

void* ReceiveRoutine(void* args){
    //work
}

int main(int argc, char **argv) {

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

    if (provided != MPI_THREAD_MULTIPLE) {
        printf("Cant make threads\n");
        MPI_Finalize();
        return 1;
    }

    tasksDistribution();

    MPI_Finalize();
}
