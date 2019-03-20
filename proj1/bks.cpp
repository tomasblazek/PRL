/*
 * algorithm: bucket sort
 * author: xblaze31
 *
 */

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>           
#include <algorithm>    /* sort */
#include <math.h>       /* log2 */

using namespace std;

#define TAG 0

bool isLess(int a, int b){
    return a < b;
}

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

int getFirstElemIndexOfTreeLevel(int level){
    int index = 0;
    for(int i = 0; i < level; i++){
        index = index + std::pow(2,i);
    }
    return index;
}

int main(int argc, char *argv[]){
    int numProcs;               //pocet procesoru
    int myid;                   //muj rank
    int neighnumber;            //hodnota souseda
    int mynumber;               //moje hodnota
    MPI_Status stat;            //struct- obsahuje kod- source, tag, error

    //MPI INIT
    MPI_Init(&argc,&argv);                          // inicializace MPI 
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);       // zjistíme, kolik procesů běží 
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);           // zjistíme id svého procesu 



    std::vector<int> vectorData;                    
    int numLeafProcs;
    int bucketSize;

    //NACTENI SOUBORU
    /* -proc s rankem 0 nacita vsechny hodnoty
     * -postupne rozesle jednotlive hodnoty vsem i sobe
    */
    if(myid == 0){
        char input[] = "numbers";                          //jmeno souboru    
        int number;                                     //hodnota pri nacitani souboru
        fstream fin;                                    //cteni ze souboru
        fin.open(input, ios::in);  

        while(fin.good()){
            number= fin.get();
            if(!fin.good()){   
                break;                      //nacte i eof, takze vyskocim
            }
            //cout<<dest<<"("<<count <<"):"<<number<<endl;             //kdo dostane kere cislo
            vectorData.push_back(number);
        }

        int countOfNumbers = vectorData.size();

        numLeafProcs = raiseToPowOf2(std::log2(countOfNumbers));
        bucketSize = countOfNumbers/numLeafProcs;
//        int maxNumbers = std::pow(2,numLeafProcs);


        if(countOfNumbers%numLeafProcs != 0){
            bucketSize++;
        }
        // cout<<numLeafProcs<<endl;
        // cout<<maxNumbers<<endl;
        // cout<<bucketSize<<endl;
        // rozeslat bucket size

        std::vector<int> metadataVector;
        metadataVector.push_back(numLeafProcs);
        metadataVector.push_back(bucketSize);

        for (int i = 0; i < numProcs; i++){
            MPI_Send(&metadataVector[0], 2, MPI_INT, i, TAG, MPI_COMM_WORLD);
        }


        // Doplneni chybejicich cisel do bucketu -1
        while(vectorData.size() < numLeafProcs * bucketSize){
            vectorData.push_back(-1);
        }

        int dest = numProcs - 1; 
        std::vector<int> vectorToSend;

        // Rozeslani cisel listovym procesorům
        for (int i = 0; i < vectorData.size(); i++){
            cout<<dest<<"("<<i<<"):"<<vectorData[i]<<endl;
            vectorToSend.push_back(vectorData[i]);
            if(((i+1)%bucketSize) == 0){    
                MPI_Send(&vectorToSend[0], bucketSize, MPI_INT, dest, TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
                vectorToSend.clear();
                dest--;
            }
        }

    }

    vectorData.resize(2);
    MPI_Recv(&vectorData[0], 2, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
    numLeafProcs = vectorData[0];
    bucketSize = vectorData[1];

    //PRIJETI HODNOTY CISLA
    //vsechny leaf procesory prijmou hodnotu a zahlasi ji



    if(numProcs - numLeafProcs <= myid && myid < numProcs){
        vectorData.resize(bucketSize);
        MPI_Recv(&vectorData[0], bucketSize, MPI_INT, 0, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
    }


    // // RAZENI
    if(numProcs - numLeafProcs <= myid && myid < numProcs){
        for (int i = 0; i < numLeafProcs; i++){
            std::sort (vectorData.begin(), vectorData.end(), isLess);
        }
        //cout <<"Proc("<<myid<<") Send data to Index parent: "<<getParentIndex(myid)<<endl;
        MPI_Send(&vectorData[0], bucketSize, MPI_INT, getParentIndex(myid), TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
    }

    // if((numProcs - numLeafProcs) <= myid && myid < numProcs){
    //     for(int i=0; i<vectorData.size(); i++){
    //         cout<<"i am:"<<myid<<" my number is:"<<vectorData[i]<<endl<<std::flush;
    //     }
    // }

    int resizer;
    int levelMax = log2(numLeafProcs);
    int level;
    int firstProcIndex;

    std::vector<int> vectorData1;
    std::vector<int> vectorData2;

    for(int j = 1; j <= levelMax; j++){
        level = levelMax - j;
        firstProcIndex = getFirstElemIndexOfTreeLevel(level); 
        if (firstProcIndex <= myid && myid < firstProcIndex + std::pow(2,level)){
            //cout<<"Recieve - procTree: "<<myid<<endl;

            resizer = std::pow(2,j-1) * bucketSize;
            vectorData1.resize(resizer);
            vectorData2.resize(resizer);


            MPI_Recv(&vectorData1[0], resizer, MPI_INT, myid*2+1, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
            MPI_Recv(&vectorData2[0], resizer, MPI_INT, myid*2+2, TAG, MPI_COMM_WORLD, &stat); //buffer,velikost,typ,rank odesilatele,tag, skupina, stat
        
            vectorData.clear();
            std::merge(vectorData1.begin(), vectorData1.end(), vectorData2.begin(), vectorData2.end(), std::back_inserter(vectorData));
            for(int i=0; i<vectorData.size(); i++){
                //cout<<"merged:"<<myid<<" my number("<< i<<") is:"<<vectorData[i]<<endl<<std::flush;
            }
            //cout<<"Send to Merge by: " << myid << " to:"<< getParentIndex(myid)<<endl<<std::flush;
            if(myid != 0){
                //cout<<"Send to Merge by: " << myid << " to:"<< getParentIndex(myid)<<endl<<std::flush;
                MPI_Send(&vectorData[0], 2*resizer, MPI_INT, getParentIndex(myid), TAG, MPI_COMM_WORLD); //buffer,velikost,typ,rank prijemce,tag,komunikacni skupina
            }
        }
    }

    // //FINALNI DISTRIBUCE VYSLEDKU-----------------------------------
    if (myid == 0){
        int eraseNegativesCount = 0; //pocet nevalidnich prvku -1 ze zacatku seznamu
        for(int i=0; i<vectorData.size(); i++){
            if (vectorData[i] < 0)
                eraseNegativesCount++;
        }

        vectorData.erase(vectorData.begin(),vectorData.begin()+eraseNegativesCount);

        for(int i=0; i<vectorData.size(); i++){
            cout<<"FINAL:"<<myid<<" my number("<< i<<") is:"<<vectorData[i]<<endl<<std::flush;
        }
    }


    MPI_Finalize(); 
    return 0;

}//main
