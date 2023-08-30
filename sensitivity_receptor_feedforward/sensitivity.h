#include <iostream>
#include <fstream>
#include<boost/array.hpp>
#include "real_value.h"

using namespace std;


int sensitivity(const vector<double> &geny, const char* fnm)
{
    int param_len = geny.size();
    const size_t rowsize = param_len;
    const size_t columnsize = 2;
    // result matrix
    vector<vector<double>> parameter_sensitivity(rowsize, vector<double>(columnsize));
    for (int i = 0; i < rowsize; ++i) {
        for (int j = 0; j < columnsize; ++j) {
            parameter_sensitivity[i][j] = 0.0;
        }
    }

    
    // analysis starts
    for(int index=0; index<param_len; ++index)
    {
        cout << index << endl;
        // scan in the direction of increasing value first
        double perturbation = 0;
        double variation = 0.5;
        while(variation > 0.001)
        {
            vector<double> new_param_set;
            new_param_set = geny;
            perturbation += variation;
            new_param_set[index] *= pow(10, perturbation);
            int validation = real_value(new_param_set, 0);
            if(!validation)
            {
                perturbation -= variation;
            }
            variation /= 2;
            cout << variation << " " << perturbation << endl;  
        }
        parameter_sensitivity[index][0] = perturbation;
        
        // scan in the direction of decreasing value then
        perturbation = 0;
        variation = 0.5;
        while(variation > 0.001)
        {
            vector<double> new_param_set;
            new_param_set = geny;
            perturbation += variation;
            new_param_set[index] *= pow(10, -perturbation);
            int validation = real_value(new_param_set, 0);
            if(!validation)
            {
                perturbation -= variation;
            }
            variation /= 2;
            cout << variation << " " << -perturbation << endl;   
        }
        parameter_sensitivity[index][1] = -perturbation;
    }


    // print results
    for (int i = 0; i < rowsize; ++i) {
        for (int j = 0; j < columnsize; ++j) 
        {
            cout << parameter_sensitivity[i][j] << " ";
        }
        cout << endl;
    }

    // print results to file
    std::ofstream myfile;
    myfile.open(fnm);
    for (int i = 0; i < rowsize; ++i) {
        for (int j = 0; j < columnsize; ++j) 
        {
            myfile << parameter_sensitivity[i][j] << " ";
        }
        myfile << endl;
    }
    myfile.close();

    return 0;   
}