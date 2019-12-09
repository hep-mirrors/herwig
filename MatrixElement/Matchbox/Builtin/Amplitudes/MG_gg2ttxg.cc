// -*- C++ -*-
//
// MG_gg2ttxg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_gg2ttxg class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15


#include "MG_gg2ttxg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 


void MG_gg2ttxg::initProc(map<string, double> & MGParams) 
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
  mME.push_back(pars->ZERO); 
  
}


void MG_gg2ttxg::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel)
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
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3], hel[4]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0]);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;    
}


void MG_gg2ttxg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate all wavefunctions
  vxxxxx(p[perm[2]], mME[0], hel[0], +1, w[0]); 
  vxxxxx(p[perm[3]], mME[1], hel[1], +1, w[1]); 
  oxxxxx(p[perm[0]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[1]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 

  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[5]); 
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[6]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->MT, pars->ZERO, w[7]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->MT, pars->ZERO, w[8]); 
  FFV1_1(w[2], w[0], pars->GC_11, pars->MT, pars->ZERO, w[9]); 
  FFV1_2(w[3], w[1], pars->GC_11, pars->MT, pars->ZERO, w[10]); 
  VVV1P0_1(w[1], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_2(w[3], w[0], pars->GC_11, pars->MT, pars->ZERO, w[12]); 
  FFV1_1(w[2], w[1], pars->GC_11, pars->MT, pars->ZERO, w[13]); 
  VVV1P0_1(w[0], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[14]); 
  VVVV1P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[15]); 
  VVVV3P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[16]); 
  VVVV4P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[17]); 

  // Calculate all amplitudes
  VVV1_0(w[5], w[6], w[4], pars->GC_10, amp[0]); 
  FFV1_0(w[3], w[7], w[5], pars->GC_11, amp[1]); 
  FFV1_0(w[8], w[2], w[5], pars->GC_11, amp[2]); 
  FFV1_0(w[10], w[9], w[4], pars->GC_11, amp[3]); 
  FFV1_0(w[3], w[9], w[11], pars->GC_11, amp[4]); 
  FFV1_0(w[8], w[9], w[1], pars->GC_11, amp[5]); 
  FFV1_0(w[12], w[13], w[4], pars->GC_11, amp[6]); 
  FFV1_0(w[12], w[2], w[11], pars->GC_11, amp[7]); 
  FFV1_0(w[12], w[7], w[1], pars->GC_11, amp[8]); 
  FFV1_0(w[3], w[13], w[14], pars->GC_11, amp[9]); 
  FFV1_0(w[10], w[2], w[14], pars->GC_11, amp[10]); 
  VVV1_0(w[14], w[1], w[6], pars->GC_10, amp[11]); 
  FFV1_0(w[8], w[13], w[0], pars->GC_11, amp[12]); 
  FFV1_0(w[10], w[7], w[0], pars->GC_11, amp[13]); 
  VVV1_0(w[0], w[11], w[6], pars->GC_10, amp[14]); 
  FFV1_0(w[3], w[2], w[15], pars->GC_11, amp[15]); 
  FFV1_0(w[3], w[2], w[16], pars->GC_11, amp[16]); 
  FFV1_0(w[3], w[2], w[17], pars->GC_11, amp[17]); 

}
