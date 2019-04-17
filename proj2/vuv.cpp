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
 * @param[in]  child  The child node id
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


bool isEdgeFoward(int edgeId, string graph){
	int countOfForwardEdges = graph.length() - 1;
	if (edgeId < countOfForwardEdges){
		return true;
	} 
	return false;
}


int getIdOfEnteringNode(int edgeId, string graph){
	int countOfNodes = graph.length();

	if (isEdgeFoward(edgeId, graph)){
		return edgeId + 1;
	} else {
		return getParentIndex(edgeId - (countOfNodes - 1) + 1);
	}
}

// int getIdOfOutputNode(int edgeId, string graph){
// 	int countOfNodes = graph.length();

// 	if (!isEdgeFoward(edgeId, graph)){
// 		return edgeId - (countOfNodes - 1) + 1;
// 	} else {
// 		return getParentIndex(edgeId + 1);
// 	}
// }


int getLeftChildOfEnteringNode(int edgeId, string graph){
	int childId = getIdOfEnteringNode(edgeId, graph) * 2 + 1;
	if (childId < graph.length()){
		return childId;
	} 
	return -1;
}

int getRightChildOfEnteringNode(int edgeId, string graph){
	int childId = getIdOfEnteringNode(edgeId, graph) * 2 + 2;
	if (childId < graph.length()){
		return childId;
	} 
	return -1;
}

int getNextEdge(int edgeId, string graph){
	int countOfNodes = graph.length();

	if(isEdgeFoward(edgeId, graph)){
		int leftChildId = getLeftChildOfEnteringNode(edgeId, graph);
		if ( leftChildId > 0){
			return (edgeId + 1) * 2;
		} else {
			return edgeId + (countOfNodes - 1);
		}
	} else {
		if (getIdOfEnteringNode(edgeId, graph) != 0 || (edgeId == (countOfNodes - 1) && countOfNodes != 2)){
			int rightChildId = getRightChildOfEnteringNode(edgeId, graph);
			if ( rightChildId > 0 && (rightChildId != (edgeId - (countOfNodes - 1) + 1))){
				return edgeId - (countOfNodes - 1) + 1;
			} else {
				return ((edgeId - (countOfNodes -1 ))/2 - 1) + (countOfNodes - 1);
			}
		}
	}

	return edgeId;
 }


int log2RoundUp(float x){
	return ceil(log2(x));
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

    if (argc != 2){
    	fprintf(stderr, "Error: Malloc error occured while creating array of succesors!\n");
    	return 1;
    }

    string input = argv[1];

    if (input.length() < 1){
    	MPI_Finalize();
    	return 0;
    }

    if (input.length() == 1){
    	cout << input[0] << ":" << 0 << endl;
    	MPI_Finalize();
    	return 0;
    }

    int countOfForwardEdges = input.length() - 1;

    // SUFFIXSUM
    // initialize succesor of edge and succesor of succesor of edge with Euler path
    int *succ = (int *) malloc(sizeof(int) * numProcs);
    if (succ == NULL){
    	fprintf(stderr, "Error: Malloc error occured while creating array of succesors!\n");
    	MPI_Finalize();
    	return 2;
    }
    succ[myid] = getNextEdge(myid, input);
    succ[succ[myid]] = getNextEdge(succ[myid], input);

    // initialize weight values of edges
  	int *val = (int *) malloc(sizeof(int) * numProcs);
  	if (val == NULL){
    	fprintf(stderr, "Error: Malloc error occured while creating array of values!\n");
    	MPI_Finalize();
    	return 2;
    }

    // initialize weights of edges
    if (succ[myid] == myid){
    	val[myid] = 0;
    } else {
    	val[myid] = isEdgeFoward(myid, input) ? -1 : 1;
    }
    // distribute weights
    MPI_Allgather(&val[myid], 1, MPI_INT, val, 1, MPI_INT, MPI_COMM_WORLD);


    //Count some of suffixe
    for (int i = 0; i < log2RoundUp(numProcs); i++){
    	val[myid] = val[myid] + val[succ[myid]];
    	succ[myid] = succ[succ[myid]];
    	MPI_Allgather(&val[myid], 1, MPI_INT, val, 1, MPI_INT, MPI_COMM_WORLD);
    	MPI_Allgather(&succ[myid], 1, MPI_INT, succ, 1, MPI_INT, MPI_COMM_WORLD);
    }

    // Correction
    val[myid] += 2;

    //Print results
    MPI_Allgather(&val[myid], 1, MPI_INT, val, 1, MPI_INT, MPI_COMM_WORLD);
    if (myid == 0){
    	cout << input[0] << ":" << 0;
    	for (int i = 1; i < input.length(); i++){
    		cout << "," << input[i] << ":" << val[i-1];
    	}
    	cout << endl;
    }
   	
   	// free resources
   	free(succ);
   	free(val);

    MPI_Finalize(); 
    return 0;

}
