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
#include <set> 

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

    /* calculate size of MIS */
    map<int, list<int> > temGraph = graph;
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
    
    int maxSize = indSet.size() - markedCover.size();

    /* find MIS */
    list<int> misTemp = indSet;
    for(int vertex: misTemp){
        auto key = graph.find(vertex);
        for (int v: key->second){
            misTemp.remove(v);
        }
    }
    indSetMax = misTemp;

    int currentSize = indSetMax.size();
    int indSize = indSet.size();
    while(currentSize < maxSize && indSize > 1){
        for(it = indSetMax.begin(); it != indSetMax.end(); ++it) {
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

            if(indSetMax.size() == maxSize){
                break;
            }
        }
        
        indSet.remove(*indSet.begin());
        currentSize = indSetMax.size();
        indSize = indSet.size();
    }

    for(it = indSetMax.begin(); it != indSetMax.end(); ++it) {
        result += to_string(*it) + " ";
    }

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
