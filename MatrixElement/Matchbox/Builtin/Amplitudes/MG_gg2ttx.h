// -*- C++ -*-
//
// MG_gg2ttx.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef MG5_Sigma_sm_gg2ttx_H
#define MG5_Sigma_sm_gg2ttx_H
//
// This is the declaration of the MG_gg2ttx class.
//

#include <complex> 
#include <vector> 
#include "Parameters_sm.h"

using namespace std; 

class MG_gg2ttx
{
  public:

    // Constructor.
    MG_gg2ttx() {}

    // Initialize process using Herwig parameters
    virtual void initProc(map<string, double> & MGParams); 

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin(vector<std::complex<double> >& amps, const vector<int>& hel); 

    // Set momenta for matrix element evaluation
    void setMomenta(vector < double * > & momenta){p = momenta;}
 
    // number of external particles 
    static const int nexternal = 4; 

    // number of colour ordered amplitudes 
    static const int namplitudes = 3; 
    std::complex<double> amp[namplitudes]; 

  private:
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 

    // number of wavefunctions used in calculate_wavefunctions
    static const int nwavefuncs = 7; 
    std::complex<double> w[nwavefuncs][18]; 

    // Pointer to the model parameters
    Parameters_sm * pars; 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta 
    vector < double * > p; 
 
}; 


#endif  // MG5_Sigma_sm_gg2ttx_H
