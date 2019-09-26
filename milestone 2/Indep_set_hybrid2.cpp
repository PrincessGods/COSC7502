#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <sstream> 
#include <iterator> 
#include <map>
#include <list> 
#include <set> 
#include <mpi.h>
#include <omp.h>

using namespace std; 

/* global variables */
const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
int myrank, mysize;

MPI_Comm nodecomm;
int nodesize, noderank;
int flag;
int *model;
MPI_Win wintable;
MPI_Aint winsize;
int windisp;

list<int> indSet;
int *localMisTemp;
int *misTemp;

map<int, list<int>> read_line(FILE* file){
    int next = 0;
    int size = 0;

    map<int, list<int>> graph;
    list<int> vertices;
    string result = ""; 

    int count = 0;
    while (next != EOF) {
        next = fgetc(file);

        if(count == 0){
            result += next;
            stringstream geek(result); 
            geek >> size;
            result = "";
        } else {
            if(next != 10 && next != 13 && next != 32){
                if(next == -1){
                    if(result == ""){
                        result = "-1";
                    }
                    stringstream geek(result);
                    int vertex = 0; 
                    geek >> vertex;
                    vertices.push_back(vertex);
                    graph.insert(pair<int, list<int>>(count - 1, vertices));
                    indSet.push_back(count - 1);
                    result = "";
                }
                result += next;
            } else if(next == 32 || next == '\n') {
                if(result == ""){
                    result = "-1";
                }

                stringstream geek(result); 
                int vertex = 0; 
                geek >> vertex;
                vertices.push_back(vertex);
                result = "";
            }
        }
        
        if(next == '\n'){
            if(count > 0){
                graph.insert(pair<int, list<int>>(count - 1, vertices));
                indSet.push_back(count - 1);
                vertices = list<int>();
            }
            count ++;
        }
    }

    return graph;
}

int initSharedMisTemp(int removeInt){
    list<int> :: iterator it; 

    MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, myrank,
              MPI_INFO_NULL, &nodecomm);
    MPI_Comm_size(nodecomm, &nodesize);
    MPI_Comm_rank(nodecomm, &noderank);

    MPI_Win_allocate_shared(indSet.size()*sizeof(int), sizeof(int),
              MPI_INFO_NULL, comm, &localMisTemp, &wintable);

    MPI_Win_get_attr(wintable, MPI_WIN_MODEL, &model, &flag);

    if (1 != flag) {
      printf("Attribute MPI_WIN_MODEL not defined\n");
    } else {
        if (MPI_WIN_UNIFIED == *model) {
            //Debug log if (myrank == 0) printf("Memory model is MPI_WIN_UNIFIED\n");
        } else {
            if (myrank == 0) printf("Memory model is *not* MPI_WIN_UNIFIED\n");
            MPI_Finalize();
        }
    }
    
    misTemp = localMisTemp;
    if (noderank != 0) {
        MPI_Win_shared_query(wintable, 0, &winsize, &windisp, &misTemp);
    }
    
    MPI_Win_fence(0, wintable);
    
    if(noderank == 0){
        int arrIndex = 0;
        for (it = indSet.begin(); it != indSet.end(); ++it) {
            misTemp[arrIndex] = *it;
            arrIndex ++;
        }
        if(removeInt != -1){
            misTemp[removeInt] = -1;
        }
    }
    
    MPI_Win_fence(0, wintable);
    return 5;
}

/* find the minimum edge cover */
int findMinCover(map<int, list<int>> graph) {
    map<int, list<int> >::iterator itr;
    list<int> :: iterator it;

    map<int, list<int> > temGraph = graph;
    int minCover = 0;
    set <int, greater <int> > markedCover;
    
    for (itr = temGraph.begin(); itr != temGraph.end(); ++itr) { 
        for(it = itr->second.begin(); it != itr->second.end(); ++it) {
            auto key = temGraph.find(*it);
            if(key->second.size() > 1 &&
                *it != -1){
                key->second.clear();
                markedCover.insert(*it);
            } else if(key->second.size() == 1 &&
                        *it != -1){
                if(markedCover.find(*it) == markedCover.end()) {
                    markedCover.insert(itr->first);
                } 
            }
        }
    }
    
    return markedCover.size();
}

/* find the maximum independent set */
void findMaxIndSet(map<int, list<int>> graph, char* input, char* output) {
    map<int, list<int>>::iterator itr;
    list<int> :: iterator it; 
    list<int> indSetMax;
    string result = "";

    MPI_Request Rrequests[mysize-1], Srequests[mysize-1];
    MPI_Status status[1];

    /* calculate size of MIS */
    int maxSize;
    if(myrank == 0){
        maxSize = indSet.size() - findMinCover(graph);
    }
    MPI_Bcast(&maxSize, mysize, MPI_INT, root, comm);
    
    /* find MIS */
    initSharedMisTemp(-1);
    int index = indSet.size()/mysize;
    int begin = myrank * index;
    int end = indSet.size();

    if(myrank != mysize - 1){
        end = begin + index;
    }

    int recvMark;
    if(myrank != 0){
        MPI_Irecv(&recvMark, 1, MPI_INT, myrank-1, myrank-1, comm, &Rrequests[myrank-1]);
        MPI_Wait(&Rrequests[myrank-1], &status[0]);
    }
    
    #pragma omp parallel
    {   
        #pragma omp for schedule(dynamic, mysize)
        for(int i = begin; i < end; i++){
            auto key = graph.find(misTemp[i]);
            if(key != graph.end()){
                for (int v:key->second){
                    misTemp[v] = -1;
                }
            }
        }
    }

    if(myrank != mysize - 1){
        recvMark = 1;
        MPI_Isend(&recvMark, 1, MPI_INT, myrank+1, myrank, comm, &Srequests[myrank]);
    }

    MPI_Barrier(comm);
    
    int indSetMaxSize = 0;
    if(myrank == 0){
        int i;
        #pragma omp parallel private(i)
        {
            #pragma omp for schedule(static) reduction (+: indSetMaxSize)
            for(i = 0; i < indSet.size(); i++){
                if(misTemp[i] != -1) {
                    indSetMaxSize++;
                }
            }
        }
    }
    MPI_Bcast(&indSetMaxSize, mysize, MPI_INT, root, comm);
    cout << "indSetMaxSize: " << indSetMaxSize << '\n';
    while(indSetMaxSize < maxSize && indSet.size() > 1){
        int arrIndex = 0;
        for(it = indSet.begin(); it != indSet.end(); ++it) {
            initSharedMisTemp(arrIndex);
            arrIndex ++;

            index = indSet.size()/mysize;
            begin = myrank * index;
            end = indSet.size();
            int removeCount = 0;
            
            if(myrank != mysize - 1){
                end = begin + index;
            }

            int recvMark;
            if(myrank != 0){
                MPI_Irecv(&recvMark, 1, MPI_INT, myrank-1, myrank-1, comm, &Rrequests[myrank-1]);
                MPI_Wait(&Rrequests[myrank-1], &status[0]);
            }

            #pragma omp parallel
            {   
                #pragma omp for schedule(dynamic, mysize)
                for(int i = begin; i < end; i++){
                    auto key = graph.find(misTemp[i]);
                    if(key != graph.end()){
                        for (int v:key->second){
                            misTemp[v] = -1;
                        }
                    }
                }
            }

            if(myrank != mysize - 1){
                recvMark = 1;
                MPI_Isend(&recvMark, 1, MPI_INT, myrank+1, myrank, comm, &Srequests[myrank]);
            }

            MPI_Barrier(comm);
            cout << "fk2: " << myrank << '\n';
            int i;
            #pragma omp parallel private(i)
            {
                #pragma omp for schedule(static) reduction (+: removeCount)
                for(i = 0; i < indSet.size(); i++){
                    if(misTemp[i] == -1) {
                        removeCount++;
                    }
                }
            }
            cout << "indSetMaxSize2: " << indSetMaxSize << '\n';
            if(indSetMaxSize < indSet.size() - removeCount){
                indSetMaxSize = indSet.size() - removeCount;
            }
            
            if(indSetMaxSize == maxSize){
                break;
            }
        }

        indSet.remove(*indSet.begin());
    }

    if(myrank == 0){
        int i;
        #pragma omp parallel private(i)
        {
            #pragma omp for schedule(static)
            for(i = 0; i < indSet.size(); i++){
                if(misTemp[i] != -1) {
                    indSetMax.push_back(misTemp[i]);
                }
            }
        }

        for(it = indSetMax.begin(); it != indSetMax.end(); ++it) {
            result += to_string(*it) + " ";
        }

        printf("The maximum independent set is: %s\n", result.c_str());
    
        //output file
        FILE* outputFile;
        outputFile = fopen(output, "w");
        
        fflush(stdout);
        fprintf(outputFile, "The maximum independent set is: %s\n", result.c_str());

        printf("Input: %s, Output: %s\n", input, output);
    }
}

int main(int argc, char** argv) {
    if (argc < 3){
        printf("Error: input arguments less then 3\n");
        return 0;
    }

    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    MPI_Comm_rank(comm, &myrank);
    MPI_Comm_size(comm, &mysize);

    //open input file
    FILE* file = fopen(argv[1], "r");

    if (file == NULL) {
        printf("Opening file failed!\n");
        return 0;
    }

    //read input file
    map<int, list<int>> input;
    input = read_line(file);

    //write output file
    findMaxIndSet(input, argv[1], argv[2]);
    
    MPI_Finalize();
    return 0;
}
