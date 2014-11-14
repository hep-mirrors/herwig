//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef HelAmps_sm_H
#define HelAmps_sm_H

#include <cmath> 
#include <complex> 

using namespace std; 

namespace MG5_sm 
{
double Sgn(double e, double f); 

void ixxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void oxxxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void sxxxxx(double p[4], int nss, std::complex<double> sc[3]); 

void vxxxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void i2xxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fi[6]);

void o2xxxx(double p[4], double fmass, int nhel, int nsf, std::complex<double>
    fo[6]);

void s2xxxx(double p[4], int nss, std::complex<double> sc[3]); 

void v2xxxx(double p[4], double vmass, int nhel, int nsv, std::complex<double>
    v[6]);

void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[]);

void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[]);

void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[]);

void VVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[]);

void VVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex);

void VVVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[]);

void VVVV4P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[]);

void VVVV3P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[]);


}  // end namespace MG5_s

#endif  // HelAmps_sm_H
