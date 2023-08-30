#include <iostream>
#include <fstream>
#include<boost/array.hpp>
#include "adaint_recovery.h"

using namespace std;

const double t1 = 5.0;
const double t2 = 10.0;
const double t3 = 15.0;

const double a1 = 3.0;
const double a2 = 5.0;
const double a3 = 10.0;

vector<double> divideByNextElement(vector<double>& input) {
    vector<double> result;
    result.reserve(input.size() - 1); 

    for (size_t i = 0; i < input.size() - 1; ++i) 
    {
        result.push_back(input[i] / input[i + 1]);
    }
    return result;
}



int real_value(const vector<double> &geny, int print)
{
    double valor;
    vector<double> periods(3);
    vector<double> amplitudes(3);
    vector<double> resultados(2);
    int resultado;
    string ff = "no_file.txt";
    const char* filename = ff.data();
    vector<vector<double>> hts = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    vector<vector<double>> rts = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    vector<double> tmp;
    vector<double> tmpf;
    vector<double> tmpff;
    periods[0] = t1;
    periods[1] = t2;
    periods[2] = t3;

    amplitudes[0] = a1;
    amplitudes[1] = a2;
    amplitudes[2] = a3;
    double min_peak_height = 0.0;
    double max_peak_height = 0.95;
    for( int i=0 ; i<periods.size() ; ++i )
        {
            for(int j=0 ; j<amplitudes.size() ; ++j)
            {
                valor=adaint_recovery(resultados, periods[i], amplitudes[j], geny, 0, filename, 1);
                
                if (valor>=50.0)
                    {
                        valor = 0.0;
                    }
                hts[i][j] = valor;
                rts[i][j] = resultados[1];    
            }
        }
    
    if (print)
    {
        for(int i=0;i<hts.size();i++){
            for(int j=0;j<hts[1].size();j++)
            {
                cout << hts[i][j] << " ";
            }
            cout << endl;
        }

        for(int i=0;i<rts.size();i++){
            for(int j=0;j<rts[1].size();j++)
            {
                cout << rts[i][j] << " ";
            }
            cout << endl;
        }
    }
    
    // Intensity scan 
    int intensity_sens = 0;
    for(int i=0;i<hts.size() ; ++i)
    {
        tmp = hts[i]; // fixed period 
        if ( std::all_of(tmp.begin(), tmp.end(), [min_peak_height](double y) { return y > min_peak_height; })  )
        {
            vector<double> result = divideByNextElement(tmp);
            if (std::all_of(result.begin(), result.end(), [max_peak_height](double y) { return y < max_peak_height; }))
            {
                intensity_sens = 1;
                break;
            }
        }
    }
    // Frequency scan 
    int frequency_sens = 0;
    if (intensity_sens)
    {
        // traspose matrix
        const size_t rows = hts.size();
        const size_t columns = hts[1].size();
        vector<vector<double>> transposed(rows, vector<double>(columns));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                transposed[j][i] = hts[i][j];
            }
        }

        vector<vector<double>> transposedrts(rows, vector<double>(columns));

        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < columns; ++j) {
                transposedrts[j][i] = rts[i][j];
            }
        }

        for(int i=0;i<transposed.size() ; ++i)
        {
            tmpf = transposed[i]; // hts fixed intensity  
            tmpff = transposedrts[i]; // rts fixed intensity 
            if ( std::all_of(tmpf.begin(), tmpf.end(), [min_peak_height](double y) { return y > min_peak_height; })  )
            {
                vector<double> resultf = divideByNextElement(tmpf);
                vector<double> resultff = divideByNextElement(tmpff);
                if ((std::all_of(resultf.begin(), resultf.end(), [max_peak_height](double y) { return y < max_peak_height; }))  && (std::all_of(resultff.begin(), resultff.end(), [max_peak_height](double y) { return y < max_peak_height; })))
                {
                    frequency_sens = 1;
                    break;
                }
            }
        }
    }
    resultado=frequency_sens*intensity_sens;
    return resultado;
}


