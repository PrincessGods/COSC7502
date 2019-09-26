#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string> 
#include <iostream> 
#include <sstream> 
#include <iterator> 
#include <map>
#include <list> 
#include <mpi.h>
#include <omp.h>

using namespace std; 

/* global variables */
const MPI_Comm comm = MPI_COMM_WORLD;
const int root = 0;
int myrank, mysize;
list<int> indSet;

map<int, list<int> > read_line(FILE* file){
    int next = 0;
    int size = 0;

    map<int, list<int> > graph;
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
                    graph.insert(pair<int, list<int> >(count - 1, vertices));
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
                graph.insert(pair<int, list<int> >(count - 1, vertices));
                indSet.push_back(count - 1);
                vertices = list<int>();
            }
            count ++;
        }
    }

    return graph;
}

/* find the minimum edge cover */
int findMinCover(map<int, list<int> > graph) {
    map<int, list<int> >::iterator itr;
    list<int> :: iterator it;

    map<int, list<int> > temGraph = graph;
    int minCover = 0;

    for (itr = temGraph.begin(); itr != temGraph.end(); ++itr) { 
        for(it = itr->second.begin(); it != itr->second.end(); ++it) {
            auto key = temGraph.find(*it);
            if(key->second.size() > 1){
                key->second.clear();
                minCover ++;
            }
        }
    }

    //cout << "rank: " << myrank << ", covers: " << minCover << '\n';
    int allsum;
    MPI_Allreduce (&minCover, &allsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return allsum;
}

/* find the maximum independent set */
void findMaxIndSet(map<int, list<int> > graph, 
                char* input, char* output) {
    map<int, list<int> >::iterator itr;
    list<int> :: iterator it; 
    list<int> indSetMax;
    string result = "";

    /* printing map */ 
    // cout << "\nThe map is : \n"; 
    // cout << "\tKEY\tELEMENT\n"; 
    // for (itr = graph.begin(); itr != graph.end(); ++itr) { 
    //     cout << '\t' << itr->first ;
    //     list <int> :: iterator it; 
    //     for(it = itr->second.begin(); it != itr->second.end(); ++it) {
    //         cout << '\t' << *it; 
    //     }
    //     cout << '\n';
    // } 
    // cout << endl;

    /* calculate size of MIS */
    int maxSize = indSet.size() - findMinCover(graph);

    /* find MIS */
    list<int> misTemp = indSet;
    for(int vertex: misTemp){
        auto key = graph.find(vertex);
        for (int v: key->second){
            misTemp.remove(v);
        }
    }
    indSetMax = misTemp;

    int currentRank = 0;
    while(indSetMax.size() < maxSize && indSet.size() > 1){
        for(it = indSet.begin(); it != indSet.end(); ++it) {
            misTemp = indSet;
            misTemp.remove(*it);

            for(int vertex: misTemp){
                auto key = graph.find(vertex);
                for (int v: key->second){
                    misTemp.remove(v);
                }
            }

            if(indSetMax.size() < misTemp.size()){
                indSetMax = misTemp;
                
            }

            int buf = 0;
            if(indSetMax.size() == maxSize){
                buf = 1;
                currentRank = myrank;
                MPI_Bcast(&buf, mysize, MPI_INT, myrank, comm);
            }

            if(buf == 1){
                break;
            }
        }

        indSet.remove(*indSet.begin());
    }

    if(myrank == currentRank){
        for(it = indSetMax.begin(); it != indSetMax.end(); ++it) {
            result += to_string(*it) + " ";
        }

        printf("From rank: %d\n", myrank);
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
    map<int, list<int> > input;
    input = read_line(file);

    //write output file
    findMaxIndSet(input, argv[1], argv[2]);
    
    MPI_Finalize();
    return 0;
}
