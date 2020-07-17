# ifndef _MULTITASKING_H_
# define _MULTITASKING_H_

using namespace std;

extern MPI_Group world_group;      // Group containing all available ranks
extern MPI_Group* task_groups;     // Array of MPI groups dedicated to the various tasks
extern MPI_Comm* task_comms;       // Array of MPI communicators dedicated to the various tasks

int define_tasking(int, int, int, int &, int &, bool);
void free_tasking(int);
    
#endif /* _MULTITASKING_H_ */
