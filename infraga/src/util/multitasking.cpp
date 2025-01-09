# ifndef _MULTITASKING_CPP_
# define _MULTITASKING_CPP_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>
#include <mpi.h>

#include "multitasking.h"

using namespace std;

MPI_Group world_group;      // Group containing all available ranks
MPI_Group* task_groups;     // Array of MPI groups dedicated to the various tasks
MPI_Comm* task_comms;       // Array of MPI communicators dedicated to the various tasks

int define_tasking(int task_cnt, int world_rank, int world_size, int & group_rank, int & group_size, bool summarize){
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    task_groups = new MPI_Group [task_cnt];
    task_comms = new MPI_Comm [task_cnt];
    
    for(int m = 0; m < task_cnt; m++){
        // For each task...
        
        // Determine how many of the world ranks will work on this task
        int rank_cnt = 0;
        for(int n = 0; n < world_size; n++){
            if(floor((n * task_cnt) / world_size) == m){
                rank_cnt++;
            }
        }
        
        // Define the rank list mapping world_rank to group_rank
        int* rank_list;
        rank_list = new int [rank_cnt];
        
        int j = 0;
        for(int n = 0; n < world_size; n++){
            if(floor((n * task_cnt) / world_size) == m){
                rank_list[j] = n;
                j++;
            }
        }
        
        // Define the task's group and communicator
        MPI_Group_incl(world_group, rank_cnt, rank_list, &task_groups[m]);
        MPI_Comm_create_group(MPI_COMM_WORLD, task_groups[m], m, &task_comms[m]);
        if (task_comms[m] != MPI_COMM_NULL) {
            MPI_Comm_rank(task_comms[m], &group_rank);
            MPI_Comm_size(task_comms[m], &group_size);
        }
        
        MPI_Barrier(MPI_COMM_WORLD);
        delete rank_list;
    }
    
    if(summarize){
        for(int n = 0; n < world_size; n++){
            if(n == world_rank){
                cout << "World rank/size: " << world_rank << "/" << world_size << '\t';
                cout << "Group rank/size: " << group_rank << "/" << group_size << '\t';
                cout << "Task ID: " << floor((world_rank * task_cnt) / world_size) << '\n';
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    
    return floor((world_rank * task_cnt) / world_size);
}

void free_tasking(int task_cnt){
    MPI_Group_free(&world_group);
    for(int m = 0; m < task_cnt; m++){
        MPI_Group_free(&task_groups[m]);
        if (task_comms[m] != MPI_COMM_NULL) {
            MPI_Comm_free(&task_comms[m]);
        }
    }
    
    delete task_groups;
    delete task_comms;
}


#endif /* _MULTITASKING_CPP_ */
