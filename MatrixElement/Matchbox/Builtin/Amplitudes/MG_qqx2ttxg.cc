// -*- C++ -*-
//
// MG_qqx2ttxg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_qqx2ttxg class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15

#include "MG_qqx2ttxg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 

void MG_qqx2ttxg::initProc(map<string, double> & MGParams) 
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

void MG_qqx2ttxg::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel, int crossed) 
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
  calculate_wavefunctions(perm, helicities[0], crossed);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;
}

void MG_qqx2ttxg::calculate_wavefunctions(const int perm[], const int hel[], int crossed)
{
  // Calculate all wavefunctions 
  // q qbar initiated
  if (crossed==1){   
    ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
    oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  }
  // q g initiated
  else if (crossed==2){
    oxxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
    oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  }
  // qbar g initiated
  else if (crossed==3){
    ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]);
    ixxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  }

  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[5]); 
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[6]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->MT, pars->ZERO, w[7]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->MT, pars->ZERO, w[8]); 
  FFV1_2(w[0], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV1_1(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[10]); 
  // Calculate all amplitudes
  VVV1_0(w[5], w[6], w[4], pars->GC_10,  amp[0]); 
  FFV1_0(w[3], w[7], w[5], pars->GC_11,  amp[1]); 
  FFV1_0(w[8], w[2], w[5], pars->GC_11,  amp[2]); 
  FFV1_0(w[9], w[1], w[6], pars->GC_11,  amp[3]); 
  FFV1_0(w[0], w[10], w[6], pars->GC_11, amp[4]); 
}
