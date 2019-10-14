//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.5.4, 2017-03-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
// and was then modified by J. Bellm.
#include "eeuuggg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm_COLOREA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > u u~ g g g WEIGHTED<=7 @1

//--------------------------------------------------------------------------
// Initialize process.

vector<int>  eeuuggg::producePermutation(double r,vector < double * > & momenta){
  static bool initialized=false;
  if (!initialized){
    initProc("param_card.dat");
    initialized=true;
  }
  setMomenta(momenta);
  sigmaKin();
  
  static const int res[6][6] = {
    {5, 6, 7, 3, 4, 0},
    {5, 7, 6, 3, 4, 0},
    {6, 5, 7, 3, 4, 0},
    {6, 7, 5, 3, 4, 0},
    {7, 5, 6, 3, 4, 0},
    {7, 6, 5, 3, 4, 0}};
  
  double jampsum=0.;
  for( int i=0;i<6;i++) jampsum+=jamp2[0][i];
  double cur=0.;
  for(int i=0;i<6;i++){
    cur+=jamp2[0][i];
    if( cur/jampsum > r )return std::vector<int>(res[i], res[i] + sizeof res[i] / sizeof res[i][0]);
  }
  //std::cout<<"producePermutation: Upps.. Something went wrong!!";
  return  std::vector<int>();
}


void eeuuggg::initProc(string param_card_name)
{
  // Instantiate the model class and set parameters that stay fixed during run
  cout<<"\nColorea: Init process eeuuggg for rearrangement (arXiv:1801.06113).";
  pars = Parameters_sm::getInstance(); 
  SLHAReader_COLOREA slha(param_card_name); 
  pars->setIndependentParameters(slha); 
  pars->setIndependentCouplings(); 
//  pars->printIndependentParameters();
//  pars->printIndependentCouplings();
  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[6]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void eeuuggg::sigmaKin()
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 
  static bool firsttime = true; 
  if (firsttime)
  {
//    pars->printDependentParameters();
//    pars->printDependentCouplings();
    firsttime = false;
  }

  // Reset color flows
  for(int i = 0; i < 6; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 128; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
//  std::complex<double> * * wfs;
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, 1, -1}, {-1, -1,
      -1, -1, -1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      1}, {-1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1}, {-1, -1, -1,
      1, -1, -1, -1}, {-1, -1, -1, 1, -1, -1, 1}, {-1, -1, -1, 1, -1, 1, -1},
      {-1, -1, -1, 1, -1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1,
      -1, 1}, {-1, -1, -1, 1, 1, 1, -1}, {-1, -1, -1, 1, 1, 1, 1}, {-1, -1, 1,
      -1, -1, -1, -1}, {-1, -1, 1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, 1, -1},
      {-1, -1, 1, -1, -1, 1, 1}, {-1, -1, 1, -1, 1, -1, -1}, {-1, -1, 1, -1, 1,
      -1, 1}, {-1, -1, 1, -1, 1, 1, -1}, {-1, -1, 1, -1, 1, 1, 1}, {-1, -1, 1,
      1, -1, -1, -1}, {-1, -1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, -1, 1, -1},
      {-1, -1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      -1, 1}, {-1, -1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1}, {-1, 1, -1,
      -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, 1}, {-1, 1, -1, -1, -1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1}, {-1, 1, -1, -1, 1,
      -1, 1}, {-1, 1, -1, -1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1}, {-1, 1, -1,
      1, -1, -1, -1}, {-1, 1, -1, 1, -1, -1, 1}, {-1, 1, -1, 1, -1, 1, -1},
      {-1, 1, -1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1, -1, -1}, {-1, 1, -1, 1, 1,
      -1, 1}, {-1, 1, -1, 1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1, 1}, {-1, 1, 1, -1,
      -1, -1, -1}, {-1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, 1, -1}, {-1,
      1, 1, -1, -1, 1, 1}, {-1, 1, 1, -1, 1, -1, -1}, {-1, 1, 1, -1, 1, -1, 1},
      {-1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1}, {-1, 1, 1, 1, -1, -1,
      -1}, {-1, 1, 1, 1, -1, -1, 1}, {-1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1,
      -1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1,
      1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1}, {1,
      -1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      1, 1}, {1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, 1, -1, 1}, {1, -1, -1,
      -1, 1, 1, -1}, {1, -1, -1, -1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1}, {1,
      -1, -1, 1, -1, -1, 1}, {1, -1, -1, 1, -1, 1, -1}, {1, -1, -1, 1, -1, 1,
      1}, {1, -1, -1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, 1,
      1, 1, -1}, {1, -1, -1, 1, 1, 1, 1}, {1, -1, 1, -1, -1, -1, -1}, {1, -1,
      1, -1, -1, -1, 1}, {1, -1, 1, -1, -1, 1, -1}, {1, -1, 1, -1, -1, 1, 1},
      {1, -1, 1, -1, 1, -1, -1}, {1, -1, 1, -1, 1, -1, 1}, {1, -1, 1, -1, 1, 1,
      -1}, {1, -1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1}, {1, -1, 1, 1,
      -1, -1, 1}, {1, -1, 1, 1, -1, 1, -1}, {1, -1, 1, 1, -1, 1, 1}, {1, -1, 1,
      1, 1, -1, -1}, {1, -1, 1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, 1, -1}, {1, -1,
      1, 1, 1, 1, 1}, {1, 1, -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, 1},
      {1, 1, -1, -1, -1, 1, -1}, {1, 1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, 1,
      -1, -1}, {1, 1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1}, {1, 1, -1, 1, -1, -1, 1}, {1, 1,
      -1, 1, -1, 1, -1}, {1, 1, -1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, -1, -1}, {1,
      1, -1, 1, 1, -1, 1}, {1, 1, -1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1}, {1,
      1, 1, -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, -1, 1, -1, -1}, {1, 1, 1, -1, 1,
      -1, 1}, {1, 1, 1, -1, 1, 1, -1}, {1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, -1,
      -1, -1}, {1, 1, 1, 1, -1, -1, 1}, {1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1,
      -1, 1, 1}, {1, 1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1,
      1, 1, -1}, {1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {24}; 

  ntry = ntry + 1; 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_1_epem_uuxggg(); 

        double tsum = 0; 
        for(int iproc = 0; iproc < nprocesses; iproc++ )
        {
          matrix_element[iproc] += t[iproc]; 
          tsum += t[iproc]; 
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel])
        {
          goodhel[ihel] = true; 
          ngood++; 
          igood[ngood] = ihel; 
        }
      }
    }
    jhel = 0; 
    sum_hel = min(sum_hel, ngood); 
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
      jhel++; 
      if (jhel >= ngood)
        jhel = 0; 
      double hwgt = double(ngood)/double(sum_hel); 
      int ihel = igood[jhel]; 
      calculate_wavefunctions(perm, helicities[ihel]); 
      t[0] = matrix_1_epem_uuxggg(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}


//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void eeuuggg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
//  int i;//, j;

  // Calculate all wavefunctions
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  FFV1P0_3(w[1], w[0], pars->GC_3, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[8]); 
  FFV1_2(w[3], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[9]); 
  FFV1_1(w[8], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[10]); 
  FFV1_1(w[8], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[11]); 
  FFV2_4_3(w[1], w[0], pars->GC_50, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[12]);
  FFV2_5_2(w[3], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[13]);
  FFV1_2(w[3], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_1(w[8], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_2(w[14], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[16]); 
  FFV2_5_1(w[8], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[17]);
  FFV2_5_2(w[14], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[18]);
  FFV1_2(w[3], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_2(w[19], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[20]); 
  FFV2_5_2(w[19], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[21]);
  VVV1P0_1(w[5], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_1(w[2], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[23]); 
  FFV1_1(w[23], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_1(w[23], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[26]); 
  FFV1_1(w[23], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[27]); 
  FFV1_2(w[26], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[28]); 
  FFV2_5_1(w[23], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[29]);
  FFV2_5_2(w[26], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[30]);
  VVV1P0_1(w[4], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_1(w[2], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_1(w[32], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[33]); 
  FFV1_1(w[32], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[34]); 
  FFV1_1(w[32], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[35]); 
  FFV2_5_1(w[32], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[36]);
  VVV1P0_1(w[4], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_1(w[2], w[7], pars->GC_2, pars->ZERO, pars->ZERO, w[38]); 
  FFV1_2(w[26], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[39]); 
  FFV1_2(w[26], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[40]); 
  FFV2_5_1(w[2], w[12], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[41]);
  FFV1_2(w[14], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  FFV1_2(w[14], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_2(w[19], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV1_2(w[19], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[45]); 
  FFV1_2(w[3], w[37], pars->GC_11, pars->ZERO, pars->ZERO, w[46]); 
  VVV1P0_1(w[37], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[47]); 
  FFV1_1(w[2], w[37], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  FFV1_2(w[3], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[49]); 
  VVV1P0_1(w[31], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[50]); 
  FFV1_1(w[2], w[31], pars->GC_11, pars->ZERO, pars->ZERO, w[51]); 
  FFV1_2(w[3], w[22], pars->GC_11, pars->ZERO, pars->ZERO, w[52]); 
  VVV1P0_1(w[4], w[22], pars->GC_10, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_1(w[2], w[22], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  VVVV1P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[55]); 
  VVVV3P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[56]); 
  VVVV4P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[57]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[9], w[10], w[6], pars->GC_11, amp[0]); 
  FFV1_0(w[9], w[11], w[5], pars->GC_11, amp[1]); 
  FFV1_0(w[13], w[10], w[6], pars->GC_11, amp[2]); 
  FFV1_0(w[13], w[11], w[5], pars->GC_11, amp[3]); 
  FFV1_0(w[14], w[15], w[6], pars->GC_11, amp[4]); 
  FFV1_0(w[16], w[8], w[6], pars->GC_11, amp[5]); 
  FFV1_0(w[14], w[17], w[6], pars->GC_11, amp[6]); 
  FFV1_0(w[18], w[8], w[6], pars->GC_11, amp[7]); 
  FFV1_0(w[19], w[15], w[5], pars->GC_11, amp[8]); 
  FFV1_0(w[20], w[8], w[5], pars->GC_11, amp[9]); 
  FFV1_0(w[19], w[17], w[5], pars->GC_11, amp[10]); 
  FFV1_0(w[21], w[8], w[5], pars->GC_11, amp[11]); 
  FFV1_0(w[3], w[15], w[22], pars->GC_11, amp[12]); 
  FFV1_0(w[9], w[8], w[22], pars->GC_11, amp[13]); 
  FFV1_0(w[3], w[17], w[22], pars->GC_11, amp[14]); 
  FFV1_0(w[13], w[8], w[22], pars->GC_11, amp[15]); 
  FFV1_0(w[9], w[24], w[6], pars->GC_11, amp[16]); 
  FFV1_0(w[9], w[25], w[4], pars->GC_11, amp[17]); 
  FFV1_0(w[13], w[24], w[6], pars->GC_11, amp[18]); 
  FFV1_0(w[13], w[25], w[4], pars->GC_11, amp[19]); 
  FFV1_0(w[26], w[27], w[6], pars->GC_11, amp[20]); 
  FFV1_0(w[28], w[23], w[6], pars->GC_11, amp[21]); 
  FFV1_0(w[26], w[29], w[6], pars->GC_11, amp[22]); 
  FFV1_0(w[30], w[23], w[6], pars->GC_11, amp[23]); 
  FFV1_0(w[19], w[27], w[4], pars->GC_11, amp[24]); 
  FFV1_0(w[20], w[23], w[4], pars->GC_11, amp[25]); 
  FFV1_0(w[19], w[29], w[4], pars->GC_11, amp[26]); 
  FFV1_0(w[21], w[23], w[4], pars->GC_11, amp[27]); 
  FFV1_0(w[3], w[27], w[31], pars->GC_11, amp[28]); 
  FFV1_0(w[9], w[23], w[31], pars->GC_11, amp[29]); 
  FFV1_0(w[3], w[29], w[31], pars->GC_11, amp[30]); 
  FFV1_0(w[13], w[23], w[31], pars->GC_11, amp[31]); 
  FFV1_0(w[9], w[33], w[5], pars->GC_11, amp[32]); 
  FFV1_0(w[9], w[34], w[4], pars->GC_11, amp[33]); 
  FFV1_0(w[13], w[33], w[5], pars->GC_11, amp[34]); 
  FFV1_0(w[13], w[34], w[4], pars->GC_11, amp[35]); 
  FFV1_0(w[26], w[35], w[5], pars->GC_11, amp[36]); 
  FFV1_0(w[28], w[32], w[5], pars->GC_11, amp[37]); 
  FFV1_0(w[26], w[36], w[5], pars->GC_11, amp[38]); 
  FFV1_0(w[30], w[32], w[5], pars->GC_11, amp[39]); 
  FFV1_0(w[14], w[35], w[4], pars->GC_11, amp[40]); 
  FFV1_0(w[16], w[32], w[4], pars->GC_11, amp[41]); 
  FFV1_0(w[14], w[36], w[4], pars->GC_11, amp[42]); 
  FFV1_0(w[18], w[32], w[4], pars->GC_11, amp[43]); 
  FFV1_0(w[3], w[35], w[37], pars->GC_11, amp[44]); 
  FFV1_0(w[9], w[32], w[37], pars->GC_11, amp[45]); 
  FFV1_0(w[3], w[36], w[37], pars->GC_11, amp[46]); 
  FFV1_0(w[13], w[32], w[37], pars->GC_11, amp[47]); 
  FFV1_0(w[39], w[38], w[6], pars->GC_11, amp[48]); 
  FFV1_0(w[40], w[38], w[5], pars->GC_11, amp[49]); 
  FFV1_0(w[39], w[41], w[6], pars->GC_11, amp[50]); 
  FFV1_0(w[40], w[41], w[5], pars->GC_11, amp[51]); 
  FFV1_0(w[26], w[38], w[22], pars->GC_11, amp[52]); 
  FFV1_0(w[28], w[2], w[22], pars->GC_11, amp[53]); 
  FFV1_0(w[26], w[41], w[22], pars->GC_11, amp[54]); 
  FFV1_0(w[30], w[2], w[22], pars->GC_11, amp[55]); 
  FFV1_0(w[42], w[38], w[6], pars->GC_11, amp[56]); 
  FFV1_0(w[43], w[38], w[4], pars->GC_11, amp[57]); 
  FFV1_0(w[42], w[41], w[6], pars->GC_11, amp[58]); 
  FFV1_0(w[43], w[41], w[4], pars->GC_11, amp[59]); 
  FFV1_0(w[14], w[38], w[31], pars->GC_11, amp[60]); 
  FFV1_0(w[16], w[2], w[31], pars->GC_11, amp[61]); 
  FFV1_0(w[14], w[41], w[31], pars->GC_11, amp[62]); 
  FFV1_0(w[18], w[2], w[31], pars->GC_11, amp[63]); 
  FFV1_0(w[44], w[38], w[5], pars->GC_11, amp[64]); 
  FFV1_0(w[45], w[38], w[4], pars->GC_11, amp[65]); 
  FFV1_0(w[44], w[41], w[5], pars->GC_11, amp[66]); 
  FFV1_0(w[45], w[41], w[4], pars->GC_11, amp[67]); 
  FFV1_0(w[19], w[38], w[37], pars->GC_11, amp[68]); 
  FFV1_0(w[20], w[2], w[37], pars->GC_11, amp[69]); 
  FFV1_0(w[19], w[41], w[37], pars->GC_11, amp[70]); 
  FFV1_0(w[21], w[2], w[37], pars->GC_11, amp[71]); 
  FFV1_0(w[46], w[38], w[6], pars->GC_11, amp[72]); 
  FFV1_0(w[3], w[38], w[47], pars->GC_11, amp[73]); 
  FFV1_0(w[9], w[48], w[6], pars->GC_11, amp[74]); 
  FFV1_0(w[9], w[2], w[47], pars->GC_11, amp[75]); 
  FFV1_0(w[46], w[41], w[6], pars->GC_11, amp[76]); 
  FFV1_0(w[3], w[41], w[47], pars->GC_11, amp[77]); 
  FFV1_0(w[13], w[48], w[6], pars->GC_11, amp[78]); 
  FFV1_0(w[13], w[2], w[47], pars->GC_11, amp[79]); 
  FFV1_0(w[49], w[38], w[5], pars->GC_11, amp[80]); 
  FFV1_0(w[3], w[38], w[50], pars->GC_11, amp[81]); 
  FFV1_0(w[9], w[51], w[5], pars->GC_11, amp[82]); 
  FFV1_0(w[9], w[2], w[50], pars->GC_11, amp[83]); 
  FFV1_0(w[49], w[41], w[5], pars->GC_11, amp[84]); 
  FFV1_0(w[3], w[41], w[50], pars->GC_11, amp[85]); 
  FFV1_0(w[13], w[51], w[5], pars->GC_11, amp[86]); 
  FFV1_0(w[13], w[2], w[50], pars->GC_11, amp[87]); 
  FFV1_0(w[52], w[38], w[4], pars->GC_11, amp[88]); 
  FFV1_0(w[3], w[38], w[53], pars->GC_11, amp[89]); 
  FFV1_0(w[9], w[54], w[4], pars->GC_11, amp[90]); 
  FFV1_0(w[9], w[2], w[53], pars->GC_11, amp[91]); 
  FFV1_0(w[52], w[41], w[4], pars->GC_11, amp[92]); 
  FFV1_0(w[3], w[41], w[53], pars->GC_11, amp[93]); 
  FFV1_0(w[13], w[54], w[4], pars->GC_11, amp[94]); 
  FFV1_0(w[13], w[2], w[53], pars->GC_11, amp[95]); 
  FFV1_0(w[3], w[38], w[55], pars->GC_11, amp[96]); 
  FFV1_0(w[3], w[38], w[56], pars->GC_11, amp[97]); 
  FFV1_0(w[3], w[38], w[57], pars->GC_11, amp[98]); 
  FFV1_0(w[9], w[2], w[55], pars->GC_11, amp[99]); 
  FFV1_0(w[9], w[2], w[56], pars->GC_11, amp[100]); 
  FFV1_0(w[9], w[2], w[57], pars->GC_11, amp[101]); 
  FFV1_0(w[3], w[41], w[55], pars->GC_11, amp[102]); 
  FFV1_0(w[3], w[41], w[56], pars->GC_11, amp[103]); 
  FFV1_0(w[3], w[41], w[57], pars->GC_11, amp[104]); 
  FFV1_0(w[13], w[2], w[55], pars->GC_11, amp[105]); 
  FFV1_0(w[13], w[2], w[56], pars->GC_11, amp[106]); 
  FFV1_0(w[13], w[2], w[57], pars->GC_11, amp[107]); 

}
double eeuuggg::matrix_1_epem_uuxggg()
{
  int i;//, j;
  // Local variables
//  const int ngraphs = 108; 
  const int ncolor = 6;
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
//  static const double denom[ncolor] = {9, 9, 9, 9, 9, 9};
//  static const double cf[ncolor][ncolor] = {{64, -8, -8, 1, 1, 10}, {-8, 64, 1,
//      10, -8, 1}, {-8, 1, 64, -8, 10, 1}, {1, 10, -8, 64, 1, -8}, {1, -8, 10,
//      1, 64, -8}, {10, 1, 1, -8, -8, 64}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[2] + amp[8] + amp[9] + amp[10] + amp[11] -
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[14] - std::complex<double>
      (0, 1) * amp[15] + amp[65] + amp[67] - std::complex<double> (0, 1) *
      amp[68] - std::complex<double> (0, 1) * amp[69] - std::complex<double>
      (0, 1) * amp[70] - std::complex<double> (0, 1) * amp[71] - amp[73] -
      std::complex<double> (0, 1) * amp[74] - amp[75] - amp[77] -
      std::complex<double> (0, 1) * amp[78] - amp[79] - std::complex<double>
      (0, 1) * amp[88] - amp[89] - amp[91] - std::complex<double> (0, 1) *
      amp[92] - amp[93] - amp[95] + amp[98] - amp[96] + amp[101] - amp[99] +
      amp[104] - amp[102] + amp[107] - amp[105];
  jamp[1] = +amp[1] + amp[3] + amp[4] + amp[5] + amp[6] + amp[7] +
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[15] + amp[57] + amp[59] - std::complex<double> (0, 1) *
      amp[60] - std::complex<double> (0, 1) * amp[61] - std::complex<double>
      (0, 1) * amp[62] - std::complex<double> (0, 1) * amp[63] - amp[81] -
      std::complex<double> (0, 1) * amp[82] - amp[83] - amp[85] -
      std::complex<double> (0, 1) * amp[86] - amp[87] + std::complex<double>
      (0, 1) * amp[88] + amp[89] + amp[91] + std::complex<double> (0, 1) *
      amp[92] + amp[93] + amp[95] + amp[96] + amp[97] + amp[99] + amp[100] +
      amp[102] + amp[103] + amp[105] + amp[106];
  jamp[2] = +amp[16] + amp[18] + amp[24] + amp[25] + amp[26] + amp[27] -
      std::complex<double> (0, 1) * amp[28] - std::complex<double> (0, 1) *
      amp[29] - std::complex<double> (0, 1) * amp[30] - std::complex<double>
      (0, 1) * amp[31] + amp[64] + amp[66] + std::complex<double> (0, 1) *
      amp[68] + std::complex<double> (0, 1) * amp[69] + std::complex<double>
      (0, 1) * amp[70] + std::complex<double> (0, 1) * amp[71] + amp[73] +
      std::complex<double> (0, 1) * amp[74] + amp[75] + amp[77] +
      std::complex<double> (0, 1) * amp[78] + amp[79] - std::complex<double>
      (0, 1) * amp[80] + amp[81] + amp[83] - std::complex<double> (0, 1) *
      amp[84] + amp[85] + amp[87] - amp[98] - amp[97] - amp[101] - amp[100] -
      amp[104] - amp[103] - amp[107] - amp[106];
  jamp[3] = +amp[17] + amp[19] + amp[20] + amp[21] + amp[22] + amp[23] +
      std::complex<double> (0, 1) * amp[28] + std::complex<double> (0, 1) *
      amp[29] + std::complex<double> (0, 1) * amp[30] + std::complex<double>
      (0, 1) * amp[31] + amp[49] + amp[51] - std::complex<double> (0, 1) *
      amp[52] - std::complex<double> (0, 1) * amp[53] - std::complex<double>
      (0, 1) * amp[54] - std::complex<double> (0, 1) * amp[55] +
      std::complex<double> (0, 1) * amp[80] - amp[81] - amp[83] +
      std::complex<double> (0, 1) * amp[84] - amp[85] - amp[87] + amp[89] -
      std::complex<double> (0, 1) * amp[90] + amp[91] + amp[93] -
      std::complex<double> (0, 1) * amp[94] + amp[95] + amp[96] + amp[97] +
      amp[99] + amp[100] + amp[102] + amp[103] + amp[105] + amp[106];
  jamp[4] = +amp[32] + amp[34] + amp[40] + amp[41] + amp[42] + amp[43] -
      std::complex<double> (0, 1) * amp[44] - std::complex<double> (0, 1) *
      amp[45] - std::complex<double> (0, 1) * amp[46] - std::complex<double>
      (0, 1) * amp[47] + amp[56] + amp[58] + std::complex<double> (0, 1) *
      amp[60] + std::complex<double> (0, 1) * amp[61] + std::complex<double>
      (0, 1) * amp[62] + std::complex<double> (0, 1) * amp[63] -
      std::complex<double> (0, 1) * amp[72] + amp[73] + amp[75] -
      std::complex<double> (0, 1) * amp[76] + amp[77] + amp[79] + amp[81] +
      std::complex<double> (0, 1) * amp[82] + amp[83] + amp[85] +
      std::complex<double> (0, 1) * amp[86] + amp[87] - amp[98] - amp[97] -
      amp[101] - amp[100] - amp[104] - amp[103] - amp[107] - amp[106];
  jamp[5] = +amp[33] + amp[35] + amp[36] + amp[37] + amp[38] + amp[39] +
      std::complex<double> (0, 1) * amp[44] + std::complex<double> (0, 1) *
      amp[45] + std::complex<double> (0, 1) * amp[46] + std::complex<double>
      (0, 1) * amp[47] + amp[48] + amp[50] + std::complex<double> (0, 1) *
      amp[52] + std::complex<double> (0, 1) * amp[53] + std::complex<double>
      (0, 1) * amp[54] + std::complex<double> (0, 1) * amp[55] +
      std::complex<double> (0, 1) * amp[72] - amp[73] - amp[75] +
      std::complex<double> (0, 1) * amp[76] - amp[77] - amp[79] - amp[89] +
      std::complex<double> (0, 1) * amp[90] - amp[91] - amp[93] +
      std::complex<double> (0, 1) * amp[94] - amp[95] + amp[98] - amp[96] +
      amp[101] - amp[99] + amp[104] - amp[102] + amp[107] - amp[105];



  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return -1.; 
}


double eeuuggg::get_jamp2(int i)
{
  return jamp2[0][i];
}

int eeuuggg::colorstring(int i, int j)
{
  static const double res[6][6] = {
    {5, 6, 7, 3, 4, 0},
    {5, 7, 6, 3, 4, 0},
    {6, 5, 7, 3, 4, 0},
    {6, 7, 5, 3, 4, 0},
    {7, 5, 6, 3, 4, 0},
    {7, 6, 5, 3, 4, 0}};
  return res[i][j];
}


int eeuuggg::NCol()
{
  const int ncolor = 6;
  return ncolor;
}





