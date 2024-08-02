#include <iostream>
#include <fstream>

#include<boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "system.h"


using namespace std;
using namespace boost::numeric::odeint;




const double ton = 1.11;

typedef boost::array< double , 6 > state_type;


double adaint(double T,  double Amax, const vector<double> &p0)
{
    vector<double> full_param;
    for( int i=0 ; i<p0.size() ; ++i )
        {
            full_param.push_back(p0[i]);
        }

    full_param.push_back(Amax);
    IFF_concat_MAX sys(full_param);
    IFF_concat_MIN sys2(full_param);


    runge_kutta4< state_type > rk4; 
    double step_size = 0.001;
    int Ton_duration = int(ton / step_size) ;
    int Toff_duration = int((T - ton)/step_size) ; 
    double max_integration_time = 5*T*10.0; 
    

    state_type x = { 0.0 , 0.0 , 0.0 , 0.0, 0.0, 0.0};
    vector<state_type> x_vec;
    vector<double> times;
    vector<double> output_variable;


    vector<double> peaks_time;
    vector<double> peaks_level; 
    int nro_picos;
    double int_threshold = 0.01;
    double t = 0.0;
    

    times.push_back(t);
    x_vec.push_back(x);
    output_variable.push_back(x[5]);
    int ht = 0;
    double min_peak_height = 0.0;
    double max_peak_height = 1.0;
    
    
    while (t <= max_integration_time)
    {
        ht+=1;
        for( size_t i=0 ; i<Ton_duration ; ++i )
        {
            rk4.do_step( sys , x , t , step_size);
            t += step_size;
            times.push_back(t);
            x_vec.push_back(x);
            output_variable.push_back(x[5]);
        }
        
        if ( (std::any_of(x.begin(), x.end(), [min_peak_height](double y) { return y < min_peak_height; })) || (std::any_of(x.begin(), x.end(), [max_peak_height](double y) { return y > max_peak_height; })) || (std::any_of(x.begin(), x.end(), [](double d) { return std::isnan(d); } )) )
        {
            return 60.0;
            break;
        }

        for( size_t i=0 ; i<Toff_duration ; ++i )
        {
            rk4.do_step( sys2 , x , t , step_size);
            t += step_size;
            times.push_back(t);
            x_vec.push_back(x);
            output_variable.push_back(x[5]);
        }

        if ( (std::any_of(x.begin(), x.end(), [min_peak_height](double y) { return y < min_peak_height; })) || (std::any_of(x.begin(), x.end(), [max_peak_height](double y) { return y > max_peak_height; })) || (std::any_of(x.begin(), x.end(), [](double d) { return std::isnan(d); } )) )
        {
            return 60.0;
            break;
        }

        // calculo maximos 
        int row = (max_element(output_variable.end()-Ton_duration-Toff_duration, output_variable.end()) - output_variable.begin());
        peaks_level.push_back(output_variable[row]);
        peaks_time.push_back(times[row]);

        nro_picos = peaks_time.size();
        if (nro_picos >= 2)
        {
            if(abs(1 - peaks_level[nro_picos-1]/peaks_level[nro_picos-2])<int_threshold)
            {
                //return (double)ht; 
                break;
            }
        }
    }
           
    return (double)ht; 
}
