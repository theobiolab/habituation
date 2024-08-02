#include <iostream>
#include <fstream>
#include<boost/array.hpp>
#include "adaint.h"

using namespace std;



double real_value(const vector<double> &geny)
{
    double valor;
    double valor_a;
    double ht_1, ht_2, ht_3;
    double ht_1_a, ht_2_a, ht_3_a;
    double diff_1, diff_2, norm;
    double diff_1_a, diff_2_a, resultado;

    ht_1 = adaint(5.0, 15.0, geny);
    ht_2 = adaint(10.0, 15.0, geny);
    ht_3 = adaint(15.0, 15.0, geny);
    if ( (ht_1>=50.0) || (ht_2>=50.0) || (ht_3>=50.0) )
        {
            valor = 0.0;
        }
    else
    {
        diff_1 = ht_1 - ht_2;
        diff_2 = ht_2 - ht_3;
        if ( (diff_1==0.0) || (diff_2==0.0)  )
        {
            norm = 2*max(abs(diff_1), abs(diff_2)) + 1.0;
            valor = (diff_1 + diff_2)/norm -1;
        }
        else
        {
            norm = 2*max(abs(diff_1), abs(diff_2));
            valor = (diff_1 + diff_2)/norm -1;
        }

    }

    ht_1_a = adaint(10.0, 10.0, geny);
    ht_2_a = adaint(10.0, 15.0, geny);
    ht_3_a = adaint(10.0, 20.0, geny);
    
    if ( (ht_1_a>=50.0) || (ht_2_a>=50.0) || (ht_3_a>=50.0) )
        {
            valor_a = 0.0;
        }
    else
    {
        diff_1_a = ht_1_a - ht_2_a;
        diff_2_a = ht_2_a - ht_3_a;
        if ( (diff_1_a==0.0) || (diff_2_a==0.0)  )
        {
            norm = 2*max(abs(diff_1_a), abs(diff_2_a)) + 1.0;
            valor_a = (diff_1_a + diff_2_a)/norm -1;
        }
        else
        {
            norm = 2*max(abs(diff_1_a), abs(diff_2_a));
            valor_a = (diff_1_a + diff_2_a)/norm -1;
        }
    }
    resultado = -abs(valor_a*valor);

    return resultado;
}
