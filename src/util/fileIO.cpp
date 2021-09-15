# ifndef _FILE_IO_CPP_
# define _FILE_IO_CPP_

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <math.h>

#include "fileIO.h"

using namespace std;

//----------------------------------------//
//----------File IO Manipulation----------//
//----------------------------------------//
int file_length(string file_name){
	ifstream file_in;	file_in.open(file_name.c_str() );
	if(!file_in.is_open()){
		return 0;
	} else {
		int count = 0;
		string line;
        
		while(!file_in.eof()) {
            getline (file_in, line);
            if(line.find("#") != 0){
                count++;
            }
		}
		file_in.close();
		return count - 1;
	}
}

int file_width(string file_name){
    ifstream file_in;
    file_in.open(file_name.c_str());
    
    if(!file_in.is_open()){
        return 0;
    } else {
        double temp;
        int count = 0;
        string line;
        
        while(!file_in.eof()) {
            getline (file_in, line);
            if(line.find("#") != 0){
                stringstream ss(line);
                while(ss >> temp){
                    count++;
                }
                break;
            }
        }
        file_in.close();
        return count;
    }
}


void file_2d_dims(string file_name, int & n1, int & n2){
    double xj, yj, fj, yprev;
    int N = file_length(file_name);
    
    ifstream file_in;    file_in.open(file_name.c_str() );
    if(file_in.is_open()){
        file_in >> xj;
        file_in >> yprev;
        file_in >> fj;
        
        for (int n = 1; n < N; n++){
            file_in >> xj;
            file_in >> yj;
            file_in >> fj;
            
            if (yj  < yprev){ n2 = n; break;}
            yprev = yj;
        }
        n1 = N / n2;
        file_in.close();
    }
}

bool string2bool (std::string v){
    return !v.empty () &&
    (strcasecmp (v.c_str (), "true") == 0 || atoi (v.c_str ()) != 0);
}

#endif /* _FILE_IO_CPP_ */
