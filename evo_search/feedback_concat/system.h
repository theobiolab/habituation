#pragma once

#include <iostream>
#include <fstream>
#include<boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;



const double Rt1 = 1.0;
const double It1 = 1.0;
const double Ot1 = 1.0;
const double Rt2 = 1.0;
const double It2 = 1.0;
const double Ot2 = 1.0;

typedef boost::array< double , 6 > state_type;


// Define the class that represents the system
class IFF_concat_MAX {
    vector<double> m_gam;
public:
    IFF_concat_MAX( vector<double> &gam ) : m_gam(gam) { }
    void operator() (const state_type &x, state_type &dxdt, const double t) const {
        dxdt[0] = m_gam[14]*m_gam[0]*(Rt1-x[0]) - m_gam[1]*x[0];
        dxdt[1] = x[2]*m_gam[2]*(It1-x[1]) - m_gam[3]*x[1];
        dxdt[2] = x[0]*m_gam[4]*(Ot1-x[2]) - x[1]*m_gam[5]*x[2]/(m_gam[6]+x[2]);
        dxdt[3] = x[2]*m_gam[12]*(Rt2-x[3]) - m_gam[13]*x[3];
        dxdt[4] = x[5]*m_gam[7]*(It2-x[4]) - m_gam[8]*x[4];
        dxdt[5] = x[3]*m_gam[9]*(Ot2-x[5]) - x[4]*m_gam[10]*x[5]/(m_gam[11]+x[5]);
    }
};

class IFF_concat_MIN {
    vector<double> m_gam;
public:
    IFF_concat_MIN( vector<double> &gam ) : m_gam(gam) { }
    void operator() (const state_type &x, state_type &dxdt, const double t) const {
        dxdt[0] =  - m_gam[1]*x[0];
        dxdt[1] = x[2]*m_gam[2]*(It1-x[1]) - m_gam[3]*x[1];
        dxdt[2] = x[0]*m_gam[4]*(Ot1-x[2]) - x[1]*m_gam[5]*x[2]/(m_gam[6]+x[2]);
        dxdt[3] = x[2]*m_gam[12]*(Rt2-x[3]) - m_gam[13]*x[3];
        dxdt[4] = x[5]*m_gam[7]*(It2-x[4]) - m_gam[8]*x[4];
        dxdt[5] = x[3]*m_gam[9]*(Ot2-x[5]) - x[4]*m_gam[10]*x[5]/(m_gam[11]+x[5]);
    }
};