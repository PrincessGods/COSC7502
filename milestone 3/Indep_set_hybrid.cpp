#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <sstream> 
#include <iterator> 
#include <map>
#include <list>

using namespace std; 

const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
int myrank, mysize;
list<int> indSet;
string result = ""; 

map<int, list<int>> read_line(FILE* file){
    int next = 0;
    int size = 0;

    map<int, list<int> > graph;
    list<int> vertices;

    int count = 0;
    while (next != EOF) {
        next = fgetc(file);

        if(count == 0){
            result += next;
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
                    int node = count - 1 + myrank;
                    while(node >= size){
                        node = node - size;
                    } 
                    indSet.push_back(node);
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
                int node = count - 1 + myrank;
                while(node >= size){
                    node = node - size;
                } 
                indSet.push_back(node);
                vertices = list<int>();
            } else {
                stringstream geek(result); 
                geek >> size;
                result = "";
            }
            count ++;
        }
    }

    return graph;
}

int findMaxIndSet(map<int, list<int>> graph) {
    MPI_Request Rrequests[mysize-1], Srequests[mysize-1];
    MPI_Status status[1];

    map<int, list<int>>::iterator itr;
    list<int> :: iterator it; 
    
    int end = indSet.size();

    /* find MIS */
    int* misTemp = (int*)(malloc(sizeof(int) * end));
    int arrIndex = 0;
    for (it = indSet.begin(); it != indSet.end(); ++it) {
        misTemp[arrIndex] = *it;
        arrIndex ++;
    }

    int indSize = 0;
    int loopCount = 0;
    
    int loopEnd = 0;
    if(end % mysize != 0){
        loopEnd = end / mysize + 1;
    } else {
        loopEnd = end / mysize;
    }

    while(loopCount != loopEnd){
        int removeCount = 0;

        for(int i = 0; i < end; i++){
            auto key = graph.find(misTemp[i]);
            if(key != graph.end()){
                for (int v:key->second){
                    if(v != -1){
                        int index = v - (mysize * loopCount) - myrank;
                        while(index < 0){
                            index = end + index;
                        }
    
                        if(misTemp[index] != -1){
                            misTemp[index] = -1;
                            removeCount ++;
                        }
                    }
                }
            }
        }
        
        if(indSize < end - removeCount){
            result = "";
            #pragma omp parallel
            {   
                int i;
                string tempResult = "";
                #pragma omp for
                for(i = 0; i < end; i++) {
                    if(misTemp[i] != -1){
                        tempResult += to_string(misTemp[i]) + " ";
                    }
                }

                #pragma omp critical
                {
                    result += tempResult;
                }
            }

            indSize = end - removeCount;
        }

        arrIndex = 0;
        for (it = indSet.begin(); it != indSet.end(); ++it) {
            int value = *it + mysize * (loopCount + 1);
            while(value >= end){
                value = value - end;
            } 
            misTemp[arrIndex] = value;
            arrIndex ++;
        }

        loopCount ++;
    }
    MPI_Barrier(comm);
    
    int othersSizes = indSize;
    if(myrank != 0){
        MPI_Irecv(&othersSizes, 1, MPI_INT, myrank-1, myrank-1, comm, &Rrequests[myrank-1]);
        MPI_Wait(&Rrequests[myrank-1], &status[0]);

        if(othersSizes < indSize){
            othersSizes = indSize;
        }
    }

    if(myrank != mysize - 1){
        MPI_Isend(&othersSizes, 1, MPI_INT, myrank+1, myrank, comm, &Srequests[myrank]);
    }

    MPI_Barrier(comm);
    if(myrank == mysize - 1){
        MPI_Isend(&othersSizes, 1, MPI_INT, 0, myrank, comm, &Srequests[myrank]);
    }

    if(myrank == 0){
        MPI_Irecv(&othersSizes, 1, MPI_INT, mysize - 1, mysize - 1, comm, &Rrequests[mysize - 1]);
        MPI_Wait(&Rrequests[mysize - 1], &status[0]);
    }

    int highestRank;
    if(othersSizes == indSize){
        highestRank = myrank;
    } else {
        highestRank = -1;
    }

    free(misTemp);
    return highestRank;
}

void writeResult(char* input, char* output){
    //printf("The maximum independent set is: %s\n", result.c_str());

    //output file
    FILE* outputFile;
    outputFile = fopen(output, "w");
    
    fflush(stdout);
    fprintf(outputFile, "The maximum independent set is: %s\n", result.c_str());

    printf("Input: %s, Output: %s\n", input, output);
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
    int resultRank = findMaxIndSet(input);
    if (myrank == resultRank){
        writeResult(argv[1], argv[2]);
    }
    
    MPI_Finalize();
    return 0;
}
