#include <iostream>
#include <fstream>
#include<boost/array.hpp>

#include "sensitivity.h"

using namespace std;



// ------------------------------------
// The main program
// ------------------------------------


int main()
{    
    vector<double> geny(14);
    ifstream testFile("/Users/Maria/Desktop/phd/habituation/sensitivity_feedback_concat/system_single.txt");    
    string line;
    
    while(getline(testFile, line)){
        std::replace(line.begin(), line.end(), ',', ' ');

        stringstream ss(line);
        ss >> geny[0];
        ss >> geny[1];
        ss >> geny[2];
        ss >> geny[3];
        ss >> geny[4];
        ss >> geny[5];
        ss >> geny[6];
        ss >> geny[7];
        ss >> geny[8];
        ss >> geny[9];
        ss >> geny[10];
        ss >> geny[11];
        ss >> geny[12];
        ss >> geny[13];

        string ff = "sensitivity_feedback_.txt";
        const char* filename = ff.data();
        int sens_analy = sensitivity(geny, filename); 
    }
    testFile.close();
    return 0;   
}



