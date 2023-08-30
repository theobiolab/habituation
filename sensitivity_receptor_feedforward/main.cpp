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
    vector<double> geny(9);
    ifstream testFile("/Users/Maria/Desktop/phd/habituation/sensitivity_receptor_feedforward/system_single.txt");    
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
        
        

        string ff = "sensitivity_receptor_feedforward_.txt";
        const char* filename = ff.data();
        int sens_analy = sensitivity(geny, filename);
    }
    testFile.close();
    return 0;   
}
