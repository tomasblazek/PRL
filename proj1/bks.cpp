/**
 * subject: PRL
 * algorithm: bucket sort
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

    std::vector<int> vectorData;                    
    int numLeafProcs;
    int bucketSize;

    // Read from file
    // Root procesor load numbers from file and send them to leaf procesors
    if(myid == 0){
        char input[] = "numbers";                       // name of file   
        int number;                                     // number value after load
        fstream fin;                                    // read from file
        fin.open(input, ios::in);  

        while(fin.good()){
            number= fin.get();
            if(!fin.good()){   
                break;                      // load eof and leave
            }

            if(vectorData.empty()){ // Print input numbers to line
                cout << number;
            } else {
                cout << " " << number;
            }
            vectorData.push_back(number);
        }
        cout << endl; // end line of input numbers

        int countOfNumbers = vectorData.size();

        // Sort for only one procesor
        if (numProcs == 1){
            t1 = MPI_Wtime();

            std::sort (vectorData.begin(), vectorData.end(), isLess);
            for(int i = 0; i < countOfNumbers; i++){
                cout << vectorData[i] << endl;
            }

            t2 = MPI_Wtime();
            if(ANALYZE_TIME){
            printf("MPI_Wtime measured time: %1.4fms\n", (t2-t1)*1000);
            }

            MPI_Finalize(); 
            return 0; // end with success
        }

        numLeafProcs = raiseToPowOf2(std::log2(countOfNumbers));
        bucketSize = countOfNumbers/numLeafProcs;


        if(countOfNumbers%numLeafProcs != 0){
            bucketSize++;
        }

        std::vector<int> metadataVector;
        metadataVector.push_back(numLeafProcs);
        metadataVector.push_back(bucketSize);

        t1 = MPI_Wtime();
        for (int i = 0; i < numProcs; i++){
            MPI_Send(&metadataVector[0], 2, MPI_INT, i, TAG, MPI_COMM_WORLD);
        }


        // Fill up buckets with -1 numbers which should be deleted at the end of sort!
        while(vectorData.size() < numLeafProcs * bucketSize){
            vectorData.push_back(-1);
        }

        int dest = numProcs - 1; 
        std::vector<int> vectorToSend;

        // Send data to leaf procesors
        for (int i = 0; i < vectorData.size(); i++){
            vectorToSend.push_back(vectorData[i]);
            if(((i+1)%bucketSize) == 0){    
                MPI_Send(&vectorToSend[0], bucketSize, MPI_INT, dest, TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
                vectorToSend.clear();
                dest--;
            }
        }

    }

    // Recieve metadata
    vectorData.resize(2);
    MPI_Recv(&vectorData[0], 2, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
    numLeafProcs = vectorData[0];
    bucketSize = vectorData[1];

    // SORT ALGORITHM
    // Sort algorithm - Sort in leafs
    if(numProcs - numLeafProcs <= myid && myid < numProcs){
        vectorData.resize(bucketSize);
        MPI_Recv(&vectorData[0], bucketSize, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat

        std::sort (vectorData.begin(), vectorData.end(), isLess);
        
        MPI_Send(&vectorData[0], bucketSize, MPI_INT, getParentIndex(myid), TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
    }
       

    int resizer;
    int levelMax = log2(numLeafProcs);
    int level;
    int firstProcIndex;

    std::vector<int> vectorData1;
    std::vector<int> vectorData2;

    // Sort algorithm - Merging
    for(int j = 1; j <= levelMax; j++){
        level = levelMax - j;
        firstProcIndex = getFirstElemIndexOfTreeLevel(level); 
        if (firstProcIndex <= myid && myid < firstProcIndex + std::pow(2,level)){
            resizer = std::pow(2,j-1) * bucketSize;
            vectorData1.resize(resizer);
            vectorData2.resize(resizer);

            MPI_Recv(&vectorData1[0], resizer, MPI_INT, myid*2+1, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
            MPI_Recv(&vectorData2[0], resizer, MPI_INT, myid*2+2, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
        
            vectorData.clear();
            std::merge(vectorData1.begin(), vectorData1.end(), vectorData2.begin(), vectorData2.end(), std::back_inserter(vectorData));

            if(myid != 0){
                MPI_Send(&vectorData[0], 2*resizer, MPI_INT, getParentIndex(myid), TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
            }
        }
    }

    t2 = MPI_Wtime();

    // Final distribution of result
    if (myid == 0){
        int eraseNegativesCount = 0; // count of nonvalid elements (-1) in vector
        for(int i=0; i<vectorData.size(); i++){
            if (vectorData[i] < 0)
                eraseNegativesCount++;
        }

        vectorData.erase(vectorData.begin(),vectorData.begin()+eraseNegativesCount);

        for(int i=0; i<vectorData.size(); i++){
            cout<<vectorData[i]<<endl;
        }

        if(ANALYZE_TIME){
        printf("MPI_Wtime measured time: %1.4fms\n", (t2-t1)*1000);
        }
    }


    MPI_Finalize(); 
    return 0;

}
