#pragma once

#include <iostream>
#include <fstream>
#include<boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;


typedef boost::array< double , 6 > state_type;


// Define the class that represents the system
class IFF_concat_MAX {
    vector<double> m_gam;
public:
    IFF_concat_MAX( vector<double> &gam ) : m_gam(gam) { }
    void operator() (const state_type &x, state_type &dxdt, const double t) const {
        dxdt[0] = m_gam[2]*(1.0-x[0]-x[5]) - m_gam[0]*m_gam[10]*(1.0-x[1]-x[5]);
        dxdt[1] = m_gam[1]*(1.0-x[0]-x[1]) - m_gam[2]*(1.0-x[0]-x[5]) + m_gam[3]*x[3]*x[5];
        dxdt[2] = x[5]*m_gam[4]*(1.0-x[2]) - m_gam[5]*x[2];
        dxdt[3] = x[4]*m_gam[8]*(1.0-x[3]) - m_gam[9]*x[3];
        dxdt[4] = x[2]*m_gam[6]*(1.0-x[4]) - m_gam[7]*x[4];
        dxdt[5] = m_gam[0]*m_gam[10]*(1.0-x[1]-x[5]) - m_gam[1]*(1.0-x[0]-x[1]) - m_gam[3]*x[3]*x[5];
    }
};

class IFF_concat_MIN {
    vector<double> m_gam;
public:
    IFF_concat_MIN( vector<double> &gam ) : m_gam(gam) { }
    void operator() (const state_type &x, state_type &dxdt, const double t) const {
        dxdt[0] = m_gam[2]*(1.0-x[0]-x[5]) ;
        dxdt[1] = m_gam[1]*(1.0-x[0]-x[1]) - m_gam[2]*(1.0-x[0]-x[5]) + m_gam[3]*x[3]*x[5];
        dxdt[2] = x[5]*m_gam[4]*(1.0-x[2]) - m_gam[5]*x[2];
        dxdt[3] = x[4]*m_gam[8]*(1.0-x[3]) - m_gam[9]*x[3];
        dxdt[4] = x[2]*m_gam[6]*(1.0-x[4]) - m_gam[7]*x[4];
        dxdt[5] =  - m_gam[1]*(1.0-x[0]-x[1]) - m_gam[3]*x[3]*x[5];
    }
};