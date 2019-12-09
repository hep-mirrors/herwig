// -*- C++ -*-
//
// MG_qqx2ttx.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_qqx2ttx class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15
 

#include "MG_qqx2ttx.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 

void MG_qqx2ttx::initProc(map<string, double> & MGParams) 
{
  // Instantiate the model class and set parameters using Herwig values
  pars = Parameters_sm::getInstance(); 
  pars->setIndependentParameters(MGParams); 
  pars->setIndependentCouplings(); 

  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MT); 
}

void MG_qqx2ttx::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel) 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset the matrix elements
  for(int i = 0; i < namplitudes; i++ )
    amp[i] = 0.; 
 
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
    perm[i] = i; 

  // Set vector of helicities
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0]);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;
}

void MG_qqx2ttx::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  swap(w[0],w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 

  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[4]); 
  // Calculate all amplitudes
  FFV1_0(w[3], w[2], w[4], pars->GC_11, amp[0]); 
}
