/**
 * subject: PRL
 * algorithm: Calculate the peak level (Výpočet úrovně vrcholu)
 * author: xblaze31
 * date: 21-03-2019
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>           
#include <algorithm>    /* sort */
#include <math.h>       /* log2 */
#include <ctime>

using namespace std;

#define TAG 0
#define ANALYZE_TIME 0 // set 1 to analyze time complexity otherwise 0

/**
 * @brief      Determines if less.
 *
 * @param[in]  a     A value
 * @param[in]  b     B value
 *
 * @return     True if less, False otherwise.
 */
bool isLess(int a, int b){
    return a < b;
}

/**
 * @brief      Gets the parent index.
 *
 * @param[in]  child  The child id
 *
 * @return     The parent index.
 */
int getParentIndex(int child){
    if (child % 2 == 1){
        return (child - 1) / 2;
    } else {
        if(child == 0)
            return -1;
        else
            return (child - 2) / 2; 
    }
}

/**
 * @brief      Raise float number to higher pow of 2.
 *
 * @param[in]  x     Float number
 *
 * @return     Integer what is pow of 2
 */
int raiseToPowOf2(double x){
    if(x < 0){
        return -1;
    }

    int i = 0;
    int pow2 = 1;
    while(x > pow2){
        i++;
        pow2 = std::pow(2,i);
    }

    return std::pow(2,i);
}

/**
 * @brief      Gets the first element index of tree level.
 *
 * @param[in]  level  The level
 *
 * @return     The first element index of tree level.
 */
int getFirstElemIndexOfTreeLevel(int level){
    int index = 0;
    for(int i = 0; i < level; i++){
        index = index + std::pow(2,i);
    }
    return index;
}

int main(int argc, char *argv[]){
    //Values for testing time 
    double t1,t2;

    int numProcs;               //count of processors
    int myid;                   //my rank
    MPI_Status stat;            //struct- source, tag, error

    //MPI INIT
    MPI_Init(&argc,&argv);                          // initialization of MPI
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);       // count of running processes
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           // get id of process




    MPI_Finalize(); 
    return 0;

}
