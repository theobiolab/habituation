#include <iostream>
#include <fstream>
#include <chrono>

#include<boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include "system.h"

//#define not_opt

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;



const double ton = 1.11;

typedef boost::array< double , 6 > state_type;

struct push_back_state_and_time
{
    std::vector< state_type >& m_states;
    std::vector< double >& m_out_states;
    std::vector< double >& m_times;

    push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &out, std::vector< double > &times )
    : m_states( states ) , m_out_states( out ) , m_times( times ) { }

    void operator()( const state_type &x_recov , double t )
    {
        m_states.push_back( x_recov );
        m_out_states.push_back( x_recov[5] );
        m_times.push_back( t );
    }
};


double adaint_recovery(vector<double> &result, double T,  double Amax, const vector<double> &p0, int print, const char* fnm, int recovery_true)
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
    double step_size = 0.01;
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
        
        integrate_const(make_controlled( 1E-12 , 1E-12 , runge_kutta_dopri5< state_type >() ) , sys , x , t , t+ton , step_size, push_back_state_and_time( x_vec ,output_variable, times ) );
        t = t+ton;
        
        integrate_const(make_controlled( 1E-12 , 1E-12 , runge_kutta_dopri5< state_type >() ) , sys2 , x , t , t+T - ton , step_size, push_back_state_and_time( x_vec ,output_variable, times ) );
        
        t = t+T - ton;

        // max element 
        int row = (max_element(output_variable.end()-Ton_duration-Toff_duration, output_variable.end()) - output_variable.begin());
        peaks_level.push_back(output_variable[row]);
        peaks_time.push_back(times[row]);

        nro_picos = peaks_time.size(); 
        if (nro_picos >= 2)
        {
            if(abs(1 - peaks_level[nro_picos-1]/peaks_level[nro_picos-2])<int_threshold)
            {
                break;
            }
        }
    }
    
    result[0] = ht-1;
    if (print)
    {
        std::ofstream myfile;
        myfile.open(fnm);
        for( size_t i=0; i<= (times.size()-1); i++ )
            {
                myfile << times[i] << " ";
                myfile << x_vec[i][0] << " ";
                myfile << x_vec[i][1] << " ";
                myfile << x_vec[i][2] << " ";
                myfile << x_vec[i][3] << " ";
                myfile << x_vec[i][4] << " ";
                myfile << x_vec[i][5] << " ";
                myfile << endl;
            }
        myfile.close();
    }
    vector<state_type> n_x_vec;
    vector<double> n_times;
    vector<double> n_output_variable;
    for( size_t i=0 ; i<times.size()-Ton_duration-Toff_duration -1; ++i )
        {
            n_times.push_back(times[i]);
            n_x_vec.push_back(x_vec[i]);
            n_output_variable.push_back(output_variable[i]);
        } 
    
    
    // ------------------------------------
    // Recovery time
    // ------------------------------------
    if ((recovery_true) && (result[0]<50))
    {
        t = n_times[n_times.size()-1];
        double tmax= T*pow(2,10) + n_times[n_times.size()-1];
        state_type x_recov = n_x_vec[n_x_vec.size()-1];
        double first_peak = peaks_level[0];

        vector<state_type> x_vec_recov; // recovery_trajectory
        vector<double> times_rec;

         
        integrate_const(make_controlled( 1E-12 , 1E-12 , runge_kutta_dopri5< state_type >() ) , sys2 , x_recov , t , tmax , step_size, push_back_state_and_time( x_vec_recov ,n_output_variable,  times_rec ) );
        

    
        // perturbation
        int dt = (int)((double)x_vec_recov.size()/2);
        int resul_t = 0;
        
        
        while (dt > 0)
        {
            
            state_type x_pert= x_vec_recov[resul_t+dt-1];
            double t_pert = 0.0;
            vector<double> output_variable_pert;
            
            integrate_const(make_controlled( 1E-12 , 1E-12 , runge_kutta_dopri5< state_type >() ) , sys , x_pert , t_pert , t_pert+ton , step_size, push_back_state_and_time( x_vec ,output_variable_pert, times ) );
            t_pert = t_pert+ton;
            
            
            integrate_const(make_controlled( 1E-12 , 1E-12 , runge_kutta_dopri5< state_type >() ) , sys2 , x_pert , t_pert , t_pert+T - ton , step_size, push_back_state_and_time( x_vec ,output_variable_pert, times ) );
        
            t_pert = t_pert+T - ton;

            // max 
            int row = (max_element(output_variable_pert.begin(), output_variable_pert.end()) - output_variable_pert.begin());
            double post_recovery_peak = output_variable_pert[row]/first_peak;
            

            if (post_recovery_peak<0.9495)
            {
                resul_t = resul_t + dt;
                
            }

            dt = (int)(dt / 2);
            
            
        }

        double recovery_time = (resul_t+1)*step_size;
        result[1] = recovery_time;
    }
    else
    {
        result[1] = -1;
    }
    return (double)(ht-1);  
}