//==========================================================================
// This file has been automatically generated for C++
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#ifndef Parameters_sm_H
#define Parameters_sm_H

#include <complex> 
#include <map>

using namespace std; 

class Parameters_sm
{
  public:

    static Parameters_sm * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double WH, WW, WZ, WT, ymtau, ymt, ymb, aS, Gf, aEWM1, MH, MZ, MTA, MT, MB,
        conjg__CKM1x1, conjg__CKM3x3, CKM3x3, MZ__exp__2, MZ__exp__4, sqrt__2,
        MH__exp__2, aEW, MW, sqrt__aEW, ee, MW__exp__2, sw__exp_2, sw2, cw, sqrt__sw2,
        sw, g1, gw, vev, vev__exp__2, lam, yb, yt, ytau, muH, ee__exp__2, sw__exp__2, 
        cw__exp__2;
    std::complex<double> complexi, I1x33, I2x33, I3x33, I4x33; 
    // Model parameters dependent on aS
    double sqrt__aS, G, G__exp__2; 
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_3, GC_50, GC_51, GC_58, GC_59; 
    // Model couplings dependent on aS
    std::complex<double> GC_10, GC_11, GC_12;

    // Set parameters that are unchanged during the run
    void setIndependentParameters(map<string, double> &MGParams); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

  private:
    static Parameters_sm * instance; 
}; 

#endif  // Parameters_sm_H

