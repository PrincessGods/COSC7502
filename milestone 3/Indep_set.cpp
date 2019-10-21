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

string read_line(FILE* file){
    
    int next = 0;
    int size = 0;

    map<int, list<int> > graph;
    list<int> vertices;
    list<int> indSet;
    list<int> indSetMax;
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

    map<int, list<int>>::iterator itr;
    list<int> :: iterator it; 

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
    
    while(loopCount != end){
        int removeCount = 0;
        for(int i = 0; i < end; i++){
            auto key = graph.find(misTemp[i]);
            if(key != graph.end()){
                for (int v:key->second){
                    if(v != -1){
                        int index = 0;
                        if(v - loopCount < 0){
                            index = end + v - loopCount;
                        } else {
                            index =  v - loopCount;
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
            for(int i = 0; i < end; i++) {
                if(misTemp[i] != -1){
                    result += to_string(misTemp[i]) + " ";
                }
            }
            indSize = end - removeCount;
        }

        arrIndex = 0;
        for (it = indSet.begin(); it != indSet.end(); ++it) {
            int value = 0;
            if(*it + loopCount + 1 >= end){
                value = *it + loopCount + 1 - end;
            } else {
                value = *it + loopCount + 1;
            }
            misTemp[arrIndex] = value;
            arrIndex ++;
        }

        loopCount ++;
    }

    free(misTemp);
    return result;
}

int main(int argc, char** argv) {
    if (argc < 3){
        printf("Error: input arguments less then 3\n");
        return 0;
    }
    //open input file
    FILE* file = fopen(argv[1], "r");

    if (file == NULL) {
        printf("Opening file failed!\n");
        return 0;
    }

    //read input file
    string input;
    input = read_line(file);
    printf("The maximum independent set is: %s\n", input.c_str());
    
    //output file
    FILE* outputFile;
    outputFile = fopen(argv[2], "w");
    
    fflush(stdout);
    fprintf(outputFile, "The maximum independent set is: %s\n", input.c_str());

    printf("Input: %s, Output: %s\n", argv[1], argv[2]);

    return 0;
}
