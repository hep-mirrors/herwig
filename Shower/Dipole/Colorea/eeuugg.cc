//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.5.4, 2017-03-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
// and was then modified by J. Bellm.
#include "eeuugg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm_COLOREA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > u u~ g g WEIGHTED<=6 @1

//--------------------------------------------------------------------------
// Initialize process.

vector<int> eeuugg::producePermutation(double r,vector < double * > & momenta){
  
  setMomenta(momenta);
  sigmaKin();
  
  static const int res[2][5] = {
    {5, 6, 3, 4, 0},
    {6, 5, 3, 4, 0}};
  
  double jampsum=0.;
  for( int i=0;i<2;i++) jampsum+=jamp2[0][i];
//  std::cout<<"\njampsum "<<jampsum<<std::flush;
  double cur=0.;

  for(int i=0;i<2;i++){
//    std::cout<<"\njamp2[0][i] "<<jamp2[0][i]<<std::flush;
    cur+=jamp2[0][i];
    if( cur/jampsum > r )return std::vector<int>(res[i], res[i] + sizeof res[i] / sizeof res[i][0]);
  }
  std::cerr<<"\nproducePermutation: Upps.. Something went wrong!!\n"<<std::flush;
  return  std::vector<int>();
}



void eeuugg::initProc(string param_card_name)
{
  cout<<"\nColorea: Init process eeuugg for rearrangement (arXiv:1801.06113).";
  // Instantiate the model class and set parameters that stay fixed during run
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
  jamp2[0] = new double[2]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void eeuugg::sigmaKin()
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
  for(int i = 0; i < 2; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 64; 
  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel;
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1},
      {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1}, {-1, -1, -1, -1, 1, 1},
      {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1}, {-1, -1, -1, 1, 1, -1},
      {-1, -1, -1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, -1, 1},
      {-1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, 1}, {-1, -1, 1, 1, -1, -1},
      {-1, -1, 1, 1, -1, 1}, {-1, -1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1}, {-1,
      1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, -1}, {-1, 1, -1, 1, -1, 1}, {-1, 1,
      -1, 1, 1, -1}, {-1, 1, -1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, 1}, {-1, 1, 1, -1, 1, -1}, {-1, 1, 1, -1, 1, 1}, {-1, 1, 1, 1,
      -1, -1}, {-1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, -1}, {-1, 1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1, -1}, {1, -1, -1, -1, -1, 1}, {1, -1, -1, -1, 1, -1},
      {1, -1, -1, -1, 1, 1}, {1, -1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1}, {1,
      -1, -1, 1, 1, -1}, {1, -1, -1, 1, 1, 1}, {1, -1, 1, -1, -1, -1}, {1, -1,
      1, -1, -1, 1}, {1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, 1}, {1, -1, 1, 1,
      -1, -1}, {1, -1, 1, 1, -1, 1}, {1, -1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, 1}, {1, 1, -1, -1, 1, -1}, {1,
      1, -1, -1, 1, 1}, {1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1}, {1, 1, -1,
      1, 1, -1}, {1, 1, -1, 1, 1, 1}, {1, 1, 1, -1, -1, -1}, {1, 1, 1, -1, -1,
      1}, {1, 1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, 1}, {1, 1, 1, 1, -1, -1}, {1,
      1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, -1}, {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {8}; 

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
        t[0] = matrix_1_epem_uuxgg(); 

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
      t[0] = matrix_1_epem_uuxgg(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 



}


//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void eeuugg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  

  // Calculate all wavefunctions
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  FFV1P0_3(w[1], w[0], pars->GC_3, pars->ZERO, pars->ZERO, w[6]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[7]); 
  FFV1_2(w[3], w[6], pars->GC_2, pars->ZERO, pars->ZERO, w[8]); 
  FFV2_4_3(w[1], w[0], pars->GC_50, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[9]);
  FFV2_5_2(w[3], w[9], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[10]);
  FFV1_2(w[3], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_1(w[2], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[2], w[6], pars->GC_2, pars->ZERO, pars->ZERO, w[14]); 
  FFV2_5_1(w[2], w[9], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[15]);
  VVV1P0_1(w[4], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[16]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[8], w[7], w[5], pars->GC_11, amp[0]); 
  FFV1_0(w[10], w[7], w[5], pars->GC_11, amp[1]); 
  FFV1_0(w[11], w[7], w[6], pars->GC_2, amp[2]); 
  FFV2_5_0(w[11], w[7], w[9], pars->GC_51, pars->GC_58, amp[3]); 
  FFV1_0(w[8], w[12], w[4], pars->GC_11, amp[4]); 
  FFV1_0(w[10], w[12], w[4], pars->GC_11, amp[5]); 
  FFV1_0(w[13], w[12], w[6], pars->GC_2, amp[6]); 
  FFV2_5_0(w[13], w[12], w[9], pars->GC_51, pars->GC_58, amp[7]); 
  FFV1_0(w[13], w[14], w[5], pars->GC_11, amp[8]); 
  FFV1_0(w[13], w[15], w[5], pars->GC_11, amp[9]); 
  FFV1_0(w[11], w[14], w[4], pars->GC_11, amp[10]); 
  FFV1_0(w[11], w[15], w[4], pars->GC_11, amp[11]); 
  FFV1_0(w[3], w[14], w[16], pars->GC_11, amp[12]); 
  FFV1_0(w[8], w[2], w[16], pars->GC_11, amp[13]); 
  FFV1_0(w[3], w[15], w[16], pars->GC_11, amp[14]); 
  FFV1_0(w[10], w[2], w[16], pars->GC_11, amp[15]); 

}
double eeuugg::matrix_1_epem_uuxgg()
{
  //int i, j;
  // Local variables
  //const int ngraphs = 16;
  //const int ncolor = 2;
  std::complex<double> jamp[2];
  // The color matrix;
  //static const double denom[ncolor] = {3, 3};
  //static const double cf[ncolor][ncolor] = {{16, -2}, {-2, 16}};

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2] + amp[3] + amp[10] + amp[11] -
      std::complex<double> (0, 1) * amp[12] - std::complex<double> (0, 1) *
      amp[13] - std::complex<double> (0, 1) * amp[14] - std::complex<double>
      (0, 1) * amp[15];
  jamp[1] = +amp[4] + amp[5] + amp[6] + amp[7] + amp[8] + amp[9] +
      std::complex<double> (0, 1) * amp[12] + std::complex<double> (0, 1) *
      amp[13] + std::complex<double> (0, 1) * amp[14] + std::complex<double>
      (0, 1) * amp[15];

  // Store the leading color flows for choice of color
  for(int i = 0; i < 2; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return -1.;
}


double eeuugg::get_jamp2(int i)
{
  return jamp2[0][i];
}

int eeuugg::colorstring(int i, int j)
{
  static const double res[2][5] = {
    {5, 6, 3, 4, 0},
    {6, 5, 3, 4, 0}};
  return res[i][j];
}


int eeuugg::NCol()
{
  const int ncolor = 2;
  return ncolor;
}






