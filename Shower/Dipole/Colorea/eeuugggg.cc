//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.5.4, 2017-03-28
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
// and was then modified by J. Bellm.
#include "eeuugggg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm_COLOREA; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: e+ e- > u u~ g g g g WEIGHTED<=8 @1

//--------------------------------------------------------------------------
// Initialize process.

vector<int>  eeuugggg::producePermutation(double r,vector < double * > & momenta){
  static bool initialized=false;
  if (!initialized){
    initProc("param_card.dat");
    initialized=true;
  }
  setMomenta(momenta);
  sigmaKin();
  
  
  static const int res[24][8] = {
    {5, 6, 7, 8, 3, 4, 0, 0},
    {5, 6, 8, 7, 3, 4, 0, 0},
    {5, 7, 6, 8, 3, 4, 0, 0},
    {5, 7, 8, 6, 3, 4, 0, 0},
    {5, 8, 6, 7, 3, 4, 0, 0},
    {5, 8, 7, 6, 3, 4, 0, 0},
    {6, 5, 7, 8, 3, 4, 0, 0},
    {6, 5, 8, 7, 3, 4, 0, 0},
    {6, 7, 5, 8, 3, 4, 0, 0},
    {6, 7, 8, 5, 3, 4, 0, 0},
    {6, 8, 5, 7, 3, 4, 0, 0},
    {6, 8, 7, 5, 3, 4, 0, 0},
    {7, 5, 6, 8, 3, 4, 0, 0},
    {7, 5, 8, 6, 3, 4, 0, 0},
    {7, 6, 5, 8, 3, 4, 0, 0},
    {7, 6, 8, 5, 3, 4, 0, 0},
    {7, 8, 5, 6, 3, 4, 0, 0},
    {7, 8, 6, 5, 3, 4, 0, 0},
    {8, 5, 6, 7, 3, 4, 0, 0},
    {8, 5, 7, 6, 3, 4, 0, 0},
    {8, 6, 5, 7, 3, 4, 0, 0},
    {8, 6, 7, 5, 3, 4, 0, 0},
    {8, 7, 5, 6, 3, 4, 0, 0},
    {8, 7, 6, 5, 3, 4, 0, 0}};
  
  double jampsum=0.;
  for( int i=0;i<24;i++) jampsum+=jamp2[0][i];
  double cur=0.;
  for(int i=0;i<24;i++){
    cur+=jamp2[0][i];
    if( cur/jampsum > r )return std::vector<int>(res[i], res[i] + sizeof res[i] / sizeof res[i][0]);
  }
  //std::cout<<"producePermutation: Upps.. Something went wrong!!";
  return  std::vector<int>();
}



void eeuugggg::initProc(string param_card_name)
{
  cout<<"\nColorea: Init process eeuugggg for rearrangement (arXiv:1801.06113).";
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
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  jamp2[0] = new double[24]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

void eeuugggg::sigmaKin()
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
  for(int i = 0; i < 24; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants
  const int ncomb = 256; 
  static bool goodhel[ncomb] = {false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
//  std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1, -1, -1, -1,
      -1}, {-1, -1, -1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, -1, -1, 1, -1},
      {-1, -1, -1, -1, -1, -1, 1, 1}, {-1, -1, -1, -1, -1, 1, -1, -1}, {-1, -1,
      -1, -1, -1, 1, -1, 1}, {-1, -1, -1, -1, -1, 1, 1, -1}, {-1, -1, -1, -1,
      -1, 1, 1, 1}, {-1, -1, -1, -1, 1, -1, -1, -1}, {-1, -1, -1, -1, 1, -1,
      -1, 1}, {-1, -1, -1, -1, 1, -1, 1, -1}, {-1, -1, -1, -1, 1, -1, 1, 1},
      {-1, -1, -1, -1, 1, 1, -1, -1}, {-1, -1, -1, -1, 1, 1, -1, 1}, {-1, -1,
      -1, -1, 1, 1, 1, -1}, {-1, -1, -1, -1, 1, 1, 1, 1}, {-1, -1, -1, 1, -1,
      -1, -1, -1}, {-1, -1, -1, 1, -1, -1, -1, 1}, {-1, -1, -1, 1, -1, -1, 1,
      -1}, {-1, -1, -1, 1, -1, -1, 1, 1}, {-1, -1, -1, 1, -1, 1, -1, -1}, {-1,
      -1, -1, 1, -1, 1, -1, 1}, {-1, -1, -1, 1, -1, 1, 1, -1}, {-1, -1, -1, 1,
      -1, 1, 1, 1}, {-1, -1, -1, 1, 1, -1, -1, -1}, {-1, -1, -1, 1, 1, -1, -1,
      1}, {-1, -1, -1, 1, 1, -1, 1, -1}, {-1, -1, -1, 1, 1, -1, 1, 1}, {-1, -1,
      -1, 1, 1, 1, -1, -1}, {-1, -1, -1, 1, 1, 1, -1, 1}, {-1, -1, -1, 1, 1, 1,
      1, -1}, {-1, -1, -1, 1, 1, 1, 1, 1}, {-1, -1, 1, -1, -1, -1, -1, -1},
      {-1, -1, 1, -1, -1, -1, -1, 1}, {-1, -1, 1, -1, -1, -1, 1, -1}, {-1, -1,
      1, -1, -1, -1, 1, 1}, {-1, -1, 1, -1, -1, 1, -1, -1}, {-1, -1, 1, -1, -1,
      1, -1, 1}, {-1, -1, 1, -1, -1, 1, 1, -1}, {-1, -1, 1, -1, -1, 1, 1, 1},
      {-1, -1, 1, -1, 1, -1, -1, -1}, {-1, -1, 1, -1, 1, -1, -1, 1}, {-1, -1,
      1, -1, 1, -1, 1, -1}, {-1, -1, 1, -1, 1, -1, 1, 1}, {-1, -1, 1, -1, 1, 1,
      -1, -1}, {-1, -1, 1, -1, 1, 1, -1, 1}, {-1, -1, 1, -1, 1, 1, 1, -1}, {-1,
      -1, 1, -1, 1, 1, 1, 1}, {-1, -1, 1, 1, -1, -1, -1, -1}, {-1, -1, 1, 1,
      -1, -1, -1, 1}, {-1, -1, 1, 1, -1, -1, 1, -1}, {-1, -1, 1, 1, -1, -1, 1,
      1}, {-1, -1, 1, 1, -1, 1, -1, -1}, {-1, -1, 1, 1, -1, 1, -1, 1}, {-1, -1,
      1, 1, -1, 1, 1, -1}, {-1, -1, 1, 1, -1, 1, 1, 1}, {-1, -1, 1, 1, 1, -1,
      -1, -1}, {-1, -1, 1, 1, 1, -1, -1, 1}, {-1, -1, 1, 1, 1, -1, 1, -1}, {-1,
      -1, 1, 1, 1, -1, 1, 1}, {-1, -1, 1, 1, 1, 1, -1, -1}, {-1, -1, 1, 1, 1,
      1, -1, 1}, {-1, -1, 1, 1, 1, 1, 1, -1}, {-1, -1, 1, 1, 1, 1, 1, 1}, {-1,
      1, -1, -1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, -1, -1, 1}, {-1, 1, -1,
      -1, -1, -1, 1, -1}, {-1, 1, -1, -1, -1, -1, 1, 1}, {-1, 1, -1, -1, -1, 1,
      -1, -1}, {-1, 1, -1, -1, -1, 1, -1, 1}, {-1, 1, -1, -1, -1, 1, 1, -1},
      {-1, 1, -1, -1, -1, 1, 1, 1}, {-1, 1, -1, -1, 1, -1, -1, -1}, {-1, 1, -1,
      -1, 1, -1, -1, 1}, {-1, 1, -1, -1, 1, -1, 1, -1}, {-1, 1, -1, -1, 1, -1,
      1, 1}, {-1, 1, -1, -1, 1, 1, -1, -1}, {-1, 1, -1, -1, 1, 1, -1, 1}, {-1,
      1, -1, -1, 1, 1, 1, -1}, {-1, 1, -1, -1, 1, 1, 1, 1}, {-1, 1, -1, 1, -1,
      -1, -1, -1}, {-1, 1, -1, 1, -1, -1, -1, 1}, {-1, 1, -1, 1, -1, -1, 1,
      -1}, {-1, 1, -1, 1, -1, -1, 1, 1}, {-1, 1, -1, 1, -1, 1, -1, -1}, {-1, 1,
      -1, 1, -1, 1, -1, 1}, {-1, 1, -1, 1, -1, 1, 1, -1}, {-1, 1, -1, 1, -1, 1,
      1, 1}, {-1, 1, -1, 1, 1, -1, -1, -1}, {-1, 1, -1, 1, 1, -1, -1, 1}, {-1,
      1, -1, 1, 1, -1, 1, -1}, {-1, 1, -1, 1, 1, -1, 1, 1}, {-1, 1, -1, 1, 1,
      1, -1, -1}, {-1, 1, -1, 1, 1, 1, -1, 1}, {-1, 1, -1, 1, 1, 1, 1, -1},
      {-1, 1, -1, 1, 1, 1, 1, 1}, {-1, 1, 1, -1, -1, -1, -1, -1}, {-1, 1, 1,
      -1, -1, -1, -1, 1}, {-1, 1, 1, -1, -1, -1, 1, -1}, {-1, 1, 1, -1, -1, -1,
      1, 1}, {-1, 1, 1, -1, -1, 1, -1, -1}, {-1, 1, 1, -1, -1, 1, -1, 1}, {-1,
      1, 1, -1, -1, 1, 1, -1}, {-1, 1, 1, -1, -1, 1, 1, 1}, {-1, 1, 1, -1, 1,
      -1, -1, -1}, {-1, 1, 1, -1, 1, -1, -1, 1}, {-1, 1, 1, -1, 1, -1, 1, -1},
      {-1, 1, 1, -1, 1, -1, 1, 1}, {-1, 1, 1, -1, 1, 1, -1, -1}, {-1, 1, 1, -1,
      1, 1, -1, 1}, {-1, 1, 1, -1, 1, 1, 1, -1}, {-1, 1, 1, -1, 1, 1, 1, 1},
      {-1, 1, 1, 1, -1, -1, -1, -1}, {-1, 1, 1, 1, -1, -1, -1, 1}, {-1, 1, 1,
      1, -1, -1, 1, -1}, {-1, 1, 1, 1, -1, -1, 1, 1}, {-1, 1, 1, 1, -1, 1, -1,
      -1}, {-1, 1, 1, 1, -1, 1, -1, 1}, {-1, 1, 1, 1, -1, 1, 1, -1}, {-1, 1, 1,
      1, -1, 1, 1, 1}, {-1, 1, 1, 1, 1, -1, -1, -1}, {-1, 1, 1, 1, 1, -1, -1,
      1}, {-1, 1, 1, 1, 1, -1, 1, -1}, {-1, 1, 1, 1, 1, -1, 1, 1}, {-1, 1, 1,
      1, 1, 1, -1, -1}, {-1, 1, 1, 1, 1, 1, -1, 1}, {-1, 1, 1, 1, 1, 1, 1, -1},
      {-1, 1, 1, 1, 1, 1, 1, 1}, {1, -1, -1, -1, -1, -1, -1, -1}, {1, -1, -1,
      -1, -1, -1, -1, 1}, {1, -1, -1, -1, -1, -1, 1, -1}, {1, -1, -1, -1, -1,
      -1, 1, 1}, {1, -1, -1, -1, -1, 1, -1, -1}, {1, -1, -1, -1, -1, 1, -1, 1},
      {1, -1, -1, -1, -1, 1, 1, -1}, {1, -1, -1, -1, -1, 1, 1, 1}, {1, -1, -1,
      -1, 1, -1, -1, -1}, {1, -1, -1, -1, 1, -1, -1, 1}, {1, -1, -1, -1, 1, -1,
      1, -1}, {1, -1, -1, -1, 1, -1, 1, 1}, {1, -1, -1, -1, 1, 1, -1, -1}, {1,
      -1, -1, -1, 1, 1, -1, 1}, {1, -1, -1, -1, 1, 1, 1, -1}, {1, -1, -1, -1,
      1, 1, 1, 1}, {1, -1, -1, 1, -1, -1, -1, -1}, {1, -1, -1, 1, -1, -1, -1,
      1}, {1, -1, -1, 1, -1, -1, 1, -1}, {1, -1, -1, 1, -1, -1, 1, 1}, {1, -1,
      -1, 1, -1, 1, -1, -1}, {1, -1, -1, 1, -1, 1, -1, 1}, {1, -1, -1, 1, -1,
      1, 1, -1}, {1, -1, -1, 1, -1, 1, 1, 1}, {1, -1, -1, 1, 1, -1, -1, -1},
      {1, -1, -1, 1, 1, -1, -1, 1}, {1, -1, -1, 1, 1, -1, 1, -1}, {1, -1, -1,
      1, 1, -1, 1, 1}, {1, -1, -1, 1, 1, 1, -1, -1}, {1, -1, -1, 1, 1, 1, -1,
      1}, {1, -1, -1, 1, 1, 1, 1, -1}, {1, -1, -1, 1, 1, 1, 1, 1}, {1, -1, 1,
      -1, -1, -1, -1, -1}, {1, -1, 1, -1, -1, -1, -1, 1}, {1, -1, 1, -1, -1,
      -1, 1, -1}, {1, -1, 1, -1, -1, -1, 1, 1}, {1, -1, 1, -1, -1, 1, -1, -1},
      {1, -1, 1, -1, -1, 1, -1, 1}, {1, -1, 1, -1, -1, 1, 1, -1}, {1, -1, 1,
      -1, -1, 1, 1, 1}, {1, -1, 1, -1, 1, -1, -1, -1}, {1, -1, 1, -1, 1, -1,
      -1, 1}, {1, -1, 1, -1, 1, -1, 1, -1}, {1, -1, 1, -1, 1, -1, 1, 1}, {1,
      -1, 1, -1, 1, 1, -1, -1}, {1, -1, 1, -1, 1, 1, -1, 1}, {1, -1, 1, -1, 1,
      1, 1, -1}, {1, -1, 1, -1, 1, 1, 1, 1}, {1, -1, 1, 1, -1, -1, -1, -1}, {1,
      -1, 1, 1, -1, -1, -1, 1}, {1, -1, 1, 1, -1, -1, 1, -1}, {1, -1, 1, 1, -1,
      -1, 1, 1}, {1, -1, 1, 1, -1, 1, -1, -1}, {1, -1, 1, 1, -1, 1, -1, 1}, {1,
      -1, 1, 1, -1, 1, 1, -1}, {1, -1, 1, 1, -1, 1, 1, 1}, {1, -1, 1, 1, 1, -1,
      -1, -1}, {1, -1, 1, 1, 1, -1, -1, 1}, {1, -1, 1, 1, 1, -1, 1, -1}, {1,
      -1, 1, 1, 1, -1, 1, 1}, {1, -1, 1, 1, 1, 1, -1, -1}, {1, -1, 1, 1, 1, 1,
      -1, 1}, {1, -1, 1, 1, 1, 1, 1, -1}, {1, -1, 1, 1, 1, 1, 1, 1}, {1, 1, -1,
      -1, -1, -1, -1, -1}, {1, 1, -1, -1, -1, -1, -1, 1}, {1, 1, -1, -1, -1,
      -1, 1, -1}, {1, 1, -1, -1, -1, -1, 1, 1}, {1, 1, -1, -1, -1, 1, -1, -1},
      {1, 1, -1, -1, -1, 1, -1, 1}, {1, 1, -1, -1, -1, 1, 1, -1}, {1, 1, -1,
      -1, -1, 1, 1, 1}, {1, 1, -1, -1, 1, -1, -1, -1}, {1, 1, -1, -1, 1, -1,
      -1, 1}, {1, 1, -1, -1, 1, -1, 1, -1}, {1, 1, -1, -1, 1, -1, 1, 1}, {1, 1,
      -1, -1, 1, 1, -1, -1}, {1, 1, -1, -1, 1, 1, -1, 1}, {1, 1, -1, -1, 1, 1,
      1, -1}, {1, 1, -1, -1, 1, 1, 1, 1}, {1, 1, -1, 1, -1, -1, -1, -1}, {1, 1,
      -1, 1, -1, -1, -1, 1}, {1, 1, -1, 1, -1, -1, 1, -1}, {1, 1, -1, 1, -1,
      -1, 1, 1}, {1, 1, -1, 1, -1, 1, -1, -1}, {1, 1, -1, 1, -1, 1, -1, 1}, {1,
      1, -1, 1, -1, 1, 1, -1}, {1, 1, -1, 1, -1, 1, 1, 1}, {1, 1, -1, 1, 1, -1,
      -1, -1}, {1, 1, -1, 1, 1, -1, -1, 1}, {1, 1, -1, 1, 1, -1, 1, -1}, {1, 1,
      -1, 1, 1, -1, 1, 1}, {1, 1, -1, 1, 1, 1, -1, -1}, {1, 1, -1, 1, 1, 1, -1,
      1}, {1, 1, -1, 1, 1, 1, 1, -1}, {1, 1, -1, 1, 1, 1, 1, 1}, {1, 1, 1, -1,
      -1, -1, -1, -1}, {1, 1, 1, -1, -1, -1, -1, 1}, {1, 1, 1, -1, -1, -1, 1,
      -1}, {1, 1, 1, -1, -1, -1, 1, 1}, {1, 1, 1, -1, -1, 1, -1, -1}, {1, 1, 1,
      -1, -1, 1, -1, 1}, {1, 1, 1, -1, -1, 1, 1, -1}, {1, 1, 1, -1, -1, 1, 1,
      1}, {1, 1, 1, -1, 1, -1, -1, -1}, {1, 1, 1, -1, 1, -1, -1, 1}, {1, 1, 1,
      -1, 1, -1, 1, -1}, {1, 1, 1, -1, 1, -1, 1, 1}, {1, 1, 1, -1, 1, 1, -1,
      -1}, {1, 1, 1, -1, 1, 1, -1, 1}, {1, 1, 1, -1, 1, 1, 1, -1}, {1, 1, 1,
      -1, 1, 1, 1, 1}, {1, 1, 1, 1, -1, -1, -1, -1}, {1, 1, 1, 1, -1, -1, -1,
      1}, {1, 1, 1, 1, -1, -1, 1, -1}, {1, 1, 1, 1, -1, -1, 1, 1}, {1, 1, 1, 1,
      -1, 1, -1, -1}, {1, 1, 1, 1, -1, 1, -1, 1}, {1, 1, 1, 1, -1, 1, 1, -1},
      {1, 1, 1, 1, -1, 1, 1, 1}, {1, 1, 1, 1, 1, -1, -1, -1}, {1, 1, 1, 1, 1,
      -1, -1, 1}, {1, 1, 1, 1, 1, -1, 1, -1}, {1, 1, 1, 1, 1, -1, 1, 1}, {1, 1,
      1, 1, 1, 1, -1, -1}, {1, 1, 1, 1, 1, 1, -1, 1}, {1, 1, 1, 1, 1, 1, 1,
      -1}, {1, 1, 1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {96}; 

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
        t[0] = matrix_1_epem_uuxgggg(); 

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
      t[0] = matrix_1_epem_uuxgggg(); 

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

void eeuugggg::calculate_wavefunctions(const int perm[], const int hel[])
{

  // Calculate all wavefunctions
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  ixxxxx(p[perm[1]], mME[1], hel[1], +1, w[1]); 
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  vxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]); 
  vxxxxx(p[perm[6]], mME[6], hel[6], +1, w[6]); 
  vxxxxx(p[perm[7]], mME[7], hel[7], +1, w[7]); 
  FFV1P0_3(w[1], w[0], pars->GC_3, pars->ZERO, pars->ZERO, w[8]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV1_2(w[3], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[10]); 
  FFV1_1(w[9], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_2(w[10], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[12]); 
  FFV1_2(w[10], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[13]); 
  FFV1_1(w[9], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[14]); 
  FFV1_2(w[10], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[15]); 
  FFV1_1(w[9], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[16]); 
  FFV2_4_3(w[1], w[0], pars->GC_50, pars->GC_59, pars->mdl_MZ, pars->mdl_WZ,
      w[17]);
  FFV2_5_2(w[3], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[18]);
  FFV1_2(w[18], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[19]); 
  FFV1_2(w[18], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[20]); 
  FFV1_2(w[18], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[21]); 
  FFV1_2(w[3], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[22]); 
  FFV1_1(w[9], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[23]); 
  FFV1_2(w[22], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[24]); 
  FFV1_2(w[22], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[25]); 
  FFV1_2(w[22], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[26]); 
  FFV2_5_1(w[9], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[27]);
  FFV2_5_2(w[22], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[28]);
  VVV1P0_1(w[6], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[29]); 
  FFV1_2(w[3], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[30]); 
  FFV1_2(w[30], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[31]); 
  FFV1_2(w[30], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[32]); 
  FFV1_2(w[30], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[33]); 
  FFV2_5_2(w[30], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[34]);
  VVV1P0_1(w[5], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[35]); 
  FFV1_2(w[3], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[36]); 
  FFV1_2(w[36], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[37]); 
  FFV1_2(w[36], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[38]); 
  FFV1_2(w[36], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[39]); 
  FFV2_5_2(w[36], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[40]);
  VVV1P0_1(w[5], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[41]); 
  FFV1_2(w[3], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[42]); 
  VVV1P0_1(w[41], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[43]); 
  FFV1_1(w[9], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[44]); 
  FFV1_2(w[3], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[45]); 
  VVV1P0_1(w[35], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[46]); 
  FFV1_1(w[9], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[47]); 
  FFV1_2(w[3], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[48]); 
  VVV1P0_1(w[5], w[29], pars->GC_10, pars->ZERO, pars->ZERO, w[49]); 
  FFV1_1(w[9], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[50]); 
  VVVV1P0_1(w[5], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[51]); 
  VVVV3P0_1(w[5], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[52]); 
  VVVV4P0_1(w[5], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[53]); 
  FFV1_1(w[2], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[54]); 
  FFV1_1(w[54], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[55]); 
  FFV1_1(w[54], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[56]); 
  FFV1_2(w[10], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[57]); 
  FFV1_1(w[54], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[58]); 
  FFV1_2(w[18], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[59]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[60]); 
  FFV1_1(w[54], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[61]); 
  FFV1_2(w[60], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[62]); 
  FFV1_2(w[60], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[63]); 
  FFV1_2(w[60], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[64]); 
  FFV2_5_1(w[54], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[65]);
  FFV2_5_2(w[60], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[66]);
  FFV1_2(w[30], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[67]); 
  VVV1P0_1(w[4], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[68]); 
  FFV1_2(w[36], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[69]); 
  VVV1P0_1(w[4], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[70]); 
  FFV1_2(w[3], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[71]); 
  VVV1P0_1(w[70], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[72]); 
  FFV1_1(w[54], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[73]); 
  FFV1_2(w[3], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[74]); 
  VVV1P0_1(w[68], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[75]); 
  FFV1_1(w[54], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[76]); 
  VVV1P0_1(w[4], w[29], pars->GC_10, pars->ZERO, pars->ZERO, w[77]); 
  FFV1_1(w[54], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[78]); 
  VVVV1P0_1(w[4], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[79]); 
  VVVV3P0_1(w[4], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[80]); 
  VVVV4P0_1(w[4], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[81]); 
  FFV1_1(w[2], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[82]); 
  FFV1_1(w[82], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[83]); 
  FFV1_1(w[82], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[84]); 
  FFV1_1(w[82], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[85]); 
  FFV1_1(w[82], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[86]); 
  FFV1_2(w[60], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[87]); 
  FFV2_5_1(w[82], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[88]);
  FFV1_2(w[22], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[89]); 
  VVV1P0_1(w[4], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[90]); 
  FFV1_2(w[3], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[91]); 
  VVV1P0_1(w[90], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[92]); 
  FFV1_1(w[82], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[93]); 
  VVV1P0_1(w[68], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[94]); 
  FFV1_1(w[82], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[95]); 
  VVV1P0_1(w[4], w[35], pars->GC_10, pars->ZERO, pars->ZERO, w[96]); 
  FFV1_1(w[82], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[97]); 
  VVVV1P0_1(w[4], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[98]); 
  VVVV3P0_1(w[4], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[99]); 
  VVVV4P0_1(w[4], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[100]); 
  FFV1_1(w[2], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[101]); 
  FFV1_1(w[101], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[102]); 
  FFV1_1(w[101], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[103]); 
  FFV1_1(w[101], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[104]); 
  FFV1_1(w[101], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[105]); 
  FFV2_5_1(w[101], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[106]);
  VVV1P0_1(w[90], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[107]); 
  FFV1_1(w[101], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[108]); 
  VVV1P0_1(w[70], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[109]); 
  FFV1_1(w[101], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[110]); 
  VVV1P0_1(w[4], w[41], pars->GC_10, pars->ZERO, pars->ZERO, w[111]); 
  FFV1_1(w[101], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[112]); 
  VVVV1P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[113]); 
  VVVV3P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[114]); 
  VVVV4P0_1(w[4], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[115]); 
  FFV1_1(w[2], w[8], pars->GC_2, pars->ZERO, pars->ZERO, w[116]); 
  FFV1_1(w[116], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[117]); 
  FFV1_1(w[116], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[118]); 
  FFV1_1(w[116], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[119]); 
  FFV2_5_1(w[2], w[17], pars->GC_51, pars->GC_58, pars->ZERO, pars->ZERO,
      w[120]);
  FFV1_1(w[120], w[6], pars->GC_11, pars->ZERO, pars->ZERO, w[121]); 
  FFV1_1(w[120], w[7], pars->GC_11, pars->ZERO, pars->ZERO, w[122]); 
  FFV1_1(w[120], w[5], pars->GC_11, pars->ZERO, pars->ZERO, w[123]); 
  FFV1_2(w[60], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[124]); 
  FFV1_1(w[2], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[125]); 
  FFV1_2(w[60], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[126]); 
  FFV1_1(w[2], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[127]); 
  FFV1_2(w[60], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[128]); 
  FFV1_1(w[2], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[129]); 
  FFV1_1(w[116], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[130]); 
  FFV1_1(w[120], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[131]); 
  FFV1_2(w[22], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[132]); 
  FFV1_1(w[2], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[133]); 
  FFV1_2(w[22], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[134]); 
  FFV1_1(w[2], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[135]); 
  FFV1_2(w[22], w[29], pars->GC_11, pars->ZERO, pars->ZERO, w[136]); 
  FFV1_2(w[30], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[137]); 
  FFV1_1(w[2], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[138]); 
  FFV1_2(w[30], w[68], pars->GC_11, pars->ZERO, pars->ZERO, w[139]); 
  FFV1_2(w[30], w[35], pars->GC_11, pars->ZERO, pars->ZERO, w[140]); 
  FFV1_2(w[36], w[90], pars->GC_11, pars->ZERO, pars->ZERO, w[141]); 
  FFV1_2(w[36], w[70], pars->GC_11, pars->ZERO, pars->ZERO, w[142]); 
  FFV1_2(w[36], w[41], pars->GC_11, pars->ZERO, pars->ZERO, w[143]); 
  FFV1P0_3(w[3], w[116], pars->GC_11, pars->ZERO, pars->ZERO, w[144]); 
  VVVV1P0_1(w[90], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[145]); 
  VVVV3P0_1(w[90], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[146]); 
  VVVV4P0_1(w[90], w[6], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[147]); 
  FFV1P0_3(w[10], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[148]); 
  FFV1P0_3(w[3], w[120], pars->GC_11, pars->ZERO, pars->ZERO, w[149]); 
  FFV1P0_3(w[18], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[150]); 
  VVV1P0_1(w[90], w[29], pars->GC_10, pars->ZERO, pars->ZERO, w[151]); 
  VVVV1P0_1(w[70], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[152]); 
  VVVV3P0_1(w[70], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[153]); 
  VVVV4P0_1(w[70], w[5], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[154]); 
  VVV1P0_1(w[70], w[35], pars->GC_10, pars->ZERO, pars->ZERO, w[155]); 
  VVVV1P0_1(w[68], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[156]); 
  VVVV3P0_1(w[68], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[157]); 
  VVVV4P0_1(w[68], w[5], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[158]); 
  VVV1P0_1(w[68], w[41], pars->GC_10, pars->ZERO, pars->ZERO, w[159]); 
  VVVV1P0_1(w[4], w[41], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[160]); 
  VVVV3P0_1(w[4], w[41], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[161]); 
  VVVV4P0_1(w[4], w[41], w[7], pars->GC_12, pars->ZERO, pars->ZERO, w[162]); 
  VVVV1P0_1(w[4], w[35], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[163]); 
  VVVV3P0_1(w[4], w[35], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[164]); 
  VVVV4P0_1(w[4], w[35], w[6], pars->GC_12, pars->ZERO, pars->ZERO, w[165]); 
  VVVV1P0_1(w[4], w[5], w[29], pars->GC_12, pars->ZERO, pars->ZERO, w[166]); 
  VVVV3P0_1(w[4], w[5], w[29], pars->GC_12, pars->ZERO, pars->ZERO, w[167]); 
  VVVV4P0_1(w[4], w[5], w[29], pars->GC_12, pars->ZERO, pars->ZERO, w[168]); 
  FFV1_2(w[3], w[113], pars->GC_11, pars->ZERO, pars->ZERO, w[169]); 
  FFV1_2(w[3], w[114], pars->GC_11, pars->ZERO, pars->ZERO, w[170]); 
  FFV1_2(w[3], w[115], pars->GC_11, pars->ZERO, pars->ZERO, w[171]); 
  VVV1P0_1(w[113], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[172]); 
  VVV1P0_1(w[114], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[173]); 
  VVV1P0_1(w[115], w[7], pars->GC_10, pars->ZERO, pars->ZERO, w[174]); 
  FFV1_1(w[2], w[113], pars->GC_11, pars->ZERO, pars->ZERO, w[175]); 
  FFV1_1(w[2], w[114], pars->GC_11, pars->ZERO, pars->ZERO, w[176]); 
  FFV1_1(w[2], w[115], pars->GC_11, pars->ZERO, pars->ZERO, w[177]); 
  FFV1_2(w[3], w[98], pars->GC_11, pars->ZERO, pars->ZERO, w[178]); 
  FFV1_2(w[3], w[99], pars->GC_11, pars->ZERO, pars->ZERO, w[179]); 
  FFV1_2(w[3], w[100], pars->GC_11, pars->ZERO, pars->ZERO, w[180]); 
  VVV1P0_1(w[98], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[181]); 
  VVV1P0_1(w[99], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[182]); 
  VVV1P0_1(w[100], w[6], pars->GC_10, pars->ZERO, pars->ZERO, w[183]); 
  FFV1_1(w[2], w[98], pars->GC_11, pars->ZERO, pars->ZERO, w[184]); 
  FFV1_1(w[2], w[99], pars->GC_11, pars->ZERO, pars->ZERO, w[185]); 
  FFV1_1(w[2], w[100], pars->GC_11, pars->ZERO, pars->ZERO, w[186]); 
  FFV1_2(w[3], w[79], pars->GC_11, pars->ZERO, pars->ZERO, w[187]); 
  FFV1_2(w[3], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[188]); 
  FFV1_2(w[3], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[189]); 
  VVV1P0_1(w[79], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[190]); 
  VVV1P0_1(w[80], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[191]); 
  VVV1P0_1(w[81], w[5], pars->GC_10, pars->ZERO, pars->ZERO, w[192]); 
  FFV1_1(w[2], w[79], pars->GC_11, pars->ZERO, pars->ZERO, w[193]); 
  FFV1_1(w[2], w[80], pars->GC_11, pars->ZERO, pars->ZERO, w[194]); 
  FFV1_1(w[2], w[81], pars->GC_11, pars->ZERO, pars->ZERO, w[195]); 
  FFV1_2(w[3], w[51], pars->GC_11, pars->ZERO, pars->ZERO, w[196]); 
  FFV1_2(w[3], w[52], pars->GC_11, pars->ZERO, pars->ZERO, w[197]); 
  FFV1_2(w[3], w[53], pars->GC_11, pars->ZERO, pars->ZERO, w[198]); 
  VVV1P0_1(w[4], w[51], pars->GC_10, pars->ZERO, pars->ZERO, w[199]); 
  VVV1P0_1(w[4], w[52], pars->GC_10, pars->ZERO, pars->ZERO, w[200]); 
  VVV1P0_1(w[4], w[53], pars->GC_10, pars->ZERO, pars->ZERO, w[201]); 
  FFV1_1(w[2], w[51], pars->GC_11, pars->ZERO, pars->ZERO, w[202]); 
  FFV1_1(w[2], w[52], pars->GC_11, pars->ZERO, pars->ZERO, w[203]); 
  FFV1_1(w[2], w[53], pars->GC_11, pars->ZERO, pars->ZERO, w[204]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[12], w[11], w[7], pars->GC_11, amp[0]); 
  FFV1_0(w[13], w[11], w[6], pars->GC_11, amp[1]); 
  FFV1_0(w[15], w[14], w[7], pars->GC_11, amp[2]); 
  FFV1_0(w[13], w[14], w[5], pars->GC_11, amp[3]); 
  FFV1_0(w[15], w[16], w[6], pars->GC_11, amp[4]); 
  FFV1_0(w[12], w[16], w[5], pars->GC_11, amp[5]); 
  FFV1_0(w[19], w[11], w[7], pars->GC_11, amp[6]); 
  FFV1_0(w[20], w[11], w[6], pars->GC_11, amp[7]); 
  FFV1_0(w[21], w[14], w[7], pars->GC_11, amp[8]); 
  FFV1_0(w[20], w[14], w[5], pars->GC_11, amp[9]); 
  FFV1_0(w[21], w[16], w[6], pars->GC_11, amp[10]); 
  FFV1_0(w[19], w[16], w[5], pars->GC_11, amp[11]); 
  FFV1_0(w[24], w[23], w[7], pars->GC_11, amp[12]); 
  FFV1_0(w[25], w[23], w[6], pars->GC_11, amp[13]); 
  FFV1_0(w[26], w[14], w[7], pars->GC_11, amp[14]); 
  FFV1_0(w[26], w[16], w[6], pars->GC_11, amp[15]); 
  FFV1_0(w[25], w[14], w[8], pars->GC_2, amp[16]); 
  FFV1_0(w[24], w[16], w[8], pars->GC_2, amp[17]); 
  FFV1_0(w[24], w[27], w[7], pars->GC_11, amp[18]); 
  FFV1_0(w[25], w[27], w[6], pars->GC_11, amp[19]); 
  FFV1_0(w[28], w[14], w[7], pars->GC_11, amp[20]); 
  FFV1_0(w[28], w[16], w[6], pars->GC_11, amp[21]); 
  FFV2_5_0(w[25], w[14], w[17], pars->GC_51, pars->GC_58, amp[22]); 
  FFV2_5_0(w[24], w[16], w[17], pars->GC_51, pars->GC_58, amp[23]); 
  FFV1_0(w[22], w[23], w[29], pars->GC_11, amp[24]); 
  FFV1_0(w[26], w[9], w[29], pars->GC_11, amp[25]); 
  FFV1_0(w[22], w[27], w[29], pars->GC_11, amp[26]); 
  FFV1_0(w[28], w[9], w[29], pars->GC_11, amp[27]); 
  FFV1_0(w[31], w[23], w[7], pars->GC_11, amp[28]); 
  FFV1_0(w[32], w[23], w[5], pars->GC_11, amp[29]); 
  FFV1_0(w[33], w[11], w[7], pars->GC_11, amp[30]); 
  FFV1_0(w[33], w[16], w[5], pars->GC_11, amp[31]); 
  FFV1_0(w[32], w[11], w[8], pars->GC_2, amp[32]); 
  FFV1_0(w[31], w[16], w[8], pars->GC_2, amp[33]); 
  FFV1_0(w[31], w[27], w[7], pars->GC_11, amp[34]); 
  FFV1_0(w[32], w[27], w[5], pars->GC_11, amp[35]); 
  FFV1_0(w[34], w[11], w[7], pars->GC_11, amp[36]); 
  FFV1_0(w[34], w[16], w[5], pars->GC_11, amp[37]); 
  FFV2_5_0(w[32], w[11], w[17], pars->GC_51, pars->GC_58, amp[38]); 
  FFV2_5_0(w[31], w[16], w[17], pars->GC_51, pars->GC_58, amp[39]); 
  FFV1_0(w[30], w[23], w[35], pars->GC_11, amp[40]); 
  FFV1_0(w[33], w[9], w[35], pars->GC_11, amp[41]); 
  FFV1_0(w[30], w[27], w[35], pars->GC_11, amp[42]); 
  FFV1_0(w[34], w[9], w[35], pars->GC_11, amp[43]); 
  FFV1_0(w[37], w[23], w[6], pars->GC_11, amp[44]); 
  FFV1_0(w[38], w[23], w[5], pars->GC_11, amp[45]); 
  FFV1_0(w[39], w[11], w[6], pars->GC_11, amp[46]); 
  FFV1_0(w[39], w[14], w[5], pars->GC_11, amp[47]); 
  FFV1_0(w[38], w[11], w[8], pars->GC_2, amp[48]); 
  FFV1_0(w[37], w[14], w[8], pars->GC_2, amp[49]); 
  FFV1_0(w[37], w[27], w[6], pars->GC_11, amp[50]); 
  FFV1_0(w[38], w[27], w[5], pars->GC_11, amp[51]); 
  FFV1_0(w[40], w[11], w[6], pars->GC_11, amp[52]); 
  FFV1_0(w[40], w[14], w[5], pars->GC_11, amp[53]); 
  FFV2_5_0(w[38], w[11], w[17], pars->GC_51, pars->GC_58, amp[54]); 
  FFV2_5_0(w[37], w[14], w[17], pars->GC_51, pars->GC_58, amp[55]); 
  FFV1_0(w[36], w[23], w[41], pars->GC_11, amp[56]); 
  FFV1_0(w[39], w[9], w[41], pars->GC_11, amp[57]); 
  FFV1_0(w[36], w[27], w[41], pars->GC_11, amp[58]); 
  FFV1_0(w[40], w[9], w[41], pars->GC_11, amp[59]); 
  FFV1_0(w[42], w[23], w[7], pars->GC_11, amp[60]); 
  FFV1_0(w[3], w[23], w[43], pars->GC_11, amp[61]); 
  FFV1_0(w[10], w[44], w[7], pars->GC_11, amp[62]); 
  FFV1_0(w[10], w[16], w[41], pars->GC_11, amp[63]); 
  FFV1_0(w[10], w[9], w[43], pars->GC_11, amp[64]); 
  FFV1_0(w[42], w[16], w[8], pars->GC_2, amp[65]); 
  FFV1_0(w[42], w[27], w[7], pars->GC_11, amp[66]); 
  FFV1_0(w[3], w[27], w[43], pars->GC_11, amp[67]); 
  FFV1_0(w[18], w[44], w[7], pars->GC_11, amp[68]); 
  FFV1_0(w[18], w[16], w[41], pars->GC_11, amp[69]); 
  FFV1_0(w[18], w[9], w[43], pars->GC_11, amp[70]); 
  FFV2_5_0(w[42], w[16], w[17], pars->GC_51, pars->GC_58, amp[71]); 
  FFV1_0(w[45], w[23], w[6], pars->GC_11, amp[72]); 
  FFV1_0(w[3], w[23], w[46], pars->GC_11, amp[73]); 
  FFV1_0(w[10], w[47], w[6], pars->GC_11, amp[74]); 
  FFV1_0(w[10], w[14], w[35], pars->GC_11, amp[75]); 
  FFV1_0(w[10], w[9], w[46], pars->GC_11, amp[76]); 
  FFV1_0(w[45], w[14], w[8], pars->GC_2, amp[77]); 
  FFV1_0(w[45], w[27], w[6], pars->GC_11, amp[78]); 
  FFV1_0(w[3], w[27], w[46], pars->GC_11, amp[79]); 
  FFV1_0(w[18], w[47], w[6], pars->GC_11, amp[80]); 
  FFV1_0(w[18], w[14], w[35], pars->GC_11, amp[81]); 
  FFV1_0(w[18], w[9], w[46], pars->GC_11, amp[82]); 
  FFV2_5_0(w[45], w[14], w[17], pars->GC_51, pars->GC_58, amp[83]); 
  FFV1_0(w[48], w[23], w[5], pars->GC_11, amp[84]); 
  FFV1_0(w[3], w[23], w[49], pars->GC_11, amp[85]); 
  FFV1_0(w[10], w[11], w[29], pars->GC_11, amp[86]); 
  FFV1_0(w[10], w[50], w[5], pars->GC_11, amp[87]); 
  FFV1_0(w[10], w[9], w[49], pars->GC_11, amp[88]); 
  FFV1_0(w[48], w[11], w[8], pars->GC_2, amp[89]); 
  FFV1_0(w[48], w[27], w[5], pars->GC_11, amp[90]); 
  FFV1_0(w[3], w[27], w[49], pars->GC_11, amp[91]); 
  FFV1_0(w[18], w[11], w[29], pars->GC_11, amp[92]); 
  FFV1_0(w[18], w[50], w[5], pars->GC_11, amp[93]); 
  FFV1_0(w[18], w[9], w[49], pars->GC_11, amp[94]); 
  FFV2_5_0(w[48], w[11], w[17], pars->GC_51, pars->GC_58, amp[95]); 
  FFV1_0(w[3], w[23], w[51], pars->GC_11, amp[96]); 
  FFV1_0(w[3], w[23], w[52], pars->GC_11, amp[97]); 
  FFV1_0(w[3], w[23], w[53], pars->GC_11, amp[98]); 
  FFV1_0(w[10], w[9], w[51], pars->GC_11, amp[99]); 
  FFV1_0(w[10], w[9], w[52], pars->GC_11, amp[100]); 
  FFV1_0(w[10], w[9], w[53], pars->GC_11, amp[101]); 
  FFV1_0(w[3], w[27], w[51], pars->GC_11, amp[102]); 
  FFV1_0(w[3], w[27], w[52], pars->GC_11, amp[103]); 
  FFV1_0(w[3], w[27], w[53], pars->GC_11, amp[104]); 
  FFV1_0(w[18], w[9], w[51], pars->GC_11, amp[105]); 
  FFV1_0(w[18], w[9], w[52], pars->GC_11, amp[106]); 
  FFV1_0(w[18], w[9], w[53], pars->GC_11, amp[107]); 
  FFV1_0(w[12], w[55], w[7], pars->GC_11, amp[108]); 
  FFV1_0(w[13], w[55], w[6], pars->GC_11, amp[109]); 
  FFV1_0(w[57], w[56], w[7], pars->GC_11, amp[110]); 
  FFV1_0(w[13], w[56], w[4], pars->GC_11, amp[111]); 
  FFV1_0(w[57], w[58], w[6], pars->GC_11, amp[112]); 
  FFV1_0(w[12], w[58], w[4], pars->GC_11, amp[113]); 
  FFV1_0(w[19], w[55], w[7], pars->GC_11, amp[114]); 
  FFV1_0(w[20], w[55], w[6], pars->GC_11, amp[115]); 
  FFV1_0(w[59], w[56], w[7], pars->GC_11, amp[116]); 
  FFV1_0(w[20], w[56], w[4], pars->GC_11, amp[117]); 
  FFV1_0(w[59], w[58], w[6], pars->GC_11, amp[118]); 
  FFV1_0(w[19], w[58], w[4], pars->GC_11, amp[119]); 
  FFV1_0(w[62], w[61], w[7], pars->GC_11, amp[120]); 
  FFV1_0(w[63], w[61], w[6], pars->GC_11, amp[121]); 
  FFV1_0(w[64], w[56], w[7], pars->GC_11, amp[122]); 
  FFV1_0(w[64], w[58], w[6], pars->GC_11, amp[123]); 
  FFV1_0(w[63], w[56], w[8], pars->GC_2, amp[124]); 
  FFV1_0(w[62], w[58], w[8], pars->GC_2, amp[125]); 
  FFV1_0(w[62], w[65], w[7], pars->GC_11, amp[126]); 
  FFV1_0(w[63], w[65], w[6], pars->GC_11, amp[127]); 
  FFV1_0(w[66], w[56], w[7], pars->GC_11, amp[128]); 
  FFV1_0(w[66], w[58], w[6], pars->GC_11, amp[129]); 
  FFV2_5_0(w[63], w[56], w[17], pars->GC_51, pars->GC_58, amp[130]); 
  FFV2_5_0(w[62], w[58], w[17], pars->GC_51, pars->GC_58, amp[131]); 
  FFV1_0(w[60], w[61], w[29], pars->GC_11, amp[132]); 
  FFV1_0(w[64], w[54], w[29], pars->GC_11, amp[133]); 
  FFV1_0(w[60], w[65], w[29], pars->GC_11, amp[134]); 
  FFV1_0(w[66], w[54], w[29], pars->GC_11, amp[135]); 
  FFV1_0(w[67], w[61], w[7], pars->GC_11, amp[136]); 
  FFV1_0(w[32], w[61], w[4], pars->GC_11, amp[137]); 
  FFV1_0(w[33], w[55], w[7], pars->GC_11, amp[138]); 
  FFV1_0(w[33], w[58], w[4], pars->GC_11, amp[139]); 
  FFV1_0(w[32], w[55], w[8], pars->GC_2, amp[140]); 
  FFV1_0(w[67], w[58], w[8], pars->GC_2, amp[141]); 
  FFV1_0(w[67], w[65], w[7], pars->GC_11, amp[142]); 
  FFV1_0(w[32], w[65], w[4], pars->GC_11, amp[143]); 
  FFV1_0(w[34], w[55], w[7], pars->GC_11, amp[144]); 
  FFV1_0(w[34], w[58], w[4], pars->GC_11, amp[145]); 
  FFV2_5_0(w[32], w[55], w[17], pars->GC_51, pars->GC_58, amp[146]); 
  FFV2_5_0(w[67], w[58], w[17], pars->GC_51, pars->GC_58, amp[147]); 
  FFV1_0(w[30], w[61], w[68], pars->GC_11, amp[148]); 
  FFV1_0(w[33], w[54], w[68], pars->GC_11, amp[149]); 
  FFV1_0(w[30], w[65], w[68], pars->GC_11, amp[150]); 
  FFV1_0(w[34], w[54], w[68], pars->GC_11, amp[151]); 
  FFV1_0(w[69], w[61], w[6], pars->GC_11, amp[152]); 
  FFV1_0(w[38], w[61], w[4], pars->GC_11, amp[153]); 
  FFV1_0(w[39], w[55], w[6], pars->GC_11, amp[154]); 
  FFV1_0(w[39], w[56], w[4], pars->GC_11, amp[155]); 
  FFV1_0(w[38], w[55], w[8], pars->GC_2, amp[156]); 
  FFV1_0(w[69], w[56], w[8], pars->GC_2, amp[157]); 
  FFV1_0(w[69], w[65], w[6], pars->GC_11, amp[158]); 
  FFV1_0(w[38], w[65], w[4], pars->GC_11, amp[159]); 
  FFV1_0(w[40], w[55], w[6], pars->GC_11, amp[160]); 
  FFV1_0(w[40], w[56], w[4], pars->GC_11, amp[161]); 
  FFV2_5_0(w[38], w[55], w[17], pars->GC_51, pars->GC_58, amp[162]); 
  FFV2_5_0(w[69], w[56], w[17], pars->GC_51, pars->GC_58, amp[163]); 
  FFV1_0(w[36], w[61], w[70], pars->GC_11, amp[164]); 
  FFV1_0(w[39], w[54], w[70], pars->GC_11, amp[165]); 
  FFV1_0(w[36], w[65], w[70], pars->GC_11, amp[166]); 
  FFV1_0(w[40], w[54], w[70], pars->GC_11, amp[167]); 
  FFV1_0(w[71], w[61], w[7], pars->GC_11, amp[168]); 
  FFV1_0(w[3], w[61], w[72], pars->GC_11, amp[169]); 
  FFV1_0(w[10], w[73], w[7], pars->GC_11, amp[170]); 
  FFV1_0(w[10], w[58], w[70], pars->GC_11, amp[171]); 
  FFV1_0(w[10], w[54], w[72], pars->GC_11, amp[172]); 
  FFV1_0(w[71], w[58], w[8], pars->GC_2, amp[173]); 
  FFV1_0(w[71], w[65], w[7], pars->GC_11, amp[174]); 
  FFV1_0(w[3], w[65], w[72], pars->GC_11, amp[175]); 
  FFV1_0(w[18], w[73], w[7], pars->GC_11, amp[176]); 
  FFV1_0(w[18], w[58], w[70], pars->GC_11, amp[177]); 
  FFV1_0(w[18], w[54], w[72], pars->GC_11, amp[178]); 
  FFV2_5_0(w[71], w[58], w[17], pars->GC_51, pars->GC_58, amp[179]); 
  FFV1_0(w[74], w[61], w[6], pars->GC_11, amp[180]); 
  FFV1_0(w[3], w[61], w[75], pars->GC_11, amp[181]); 
  FFV1_0(w[10], w[76], w[6], pars->GC_11, amp[182]); 
  FFV1_0(w[10], w[56], w[68], pars->GC_11, amp[183]); 
  FFV1_0(w[10], w[54], w[75], pars->GC_11, amp[184]); 
  FFV1_0(w[74], w[56], w[8], pars->GC_2, amp[185]); 
  FFV1_0(w[74], w[65], w[6], pars->GC_11, amp[186]); 
  FFV1_0(w[3], w[65], w[75], pars->GC_11, amp[187]); 
  FFV1_0(w[18], w[76], w[6], pars->GC_11, amp[188]); 
  FFV1_0(w[18], w[56], w[68], pars->GC_11, amp[189]); 
  FFV1_0(w[18], w[54], w[75], pars->GC_11, amp[190]); 
  FFV2_5_0(w[74], w[56], w[17], pars->GC_51, pars->GC_58, amp[191]); 
  FFV1_0(w[48], w[61], w[4], pars->GC_11, amp[192]); 
  FFV1_0(w[3], w[61], w[77], pars->GC_11, amp[193]); 
  FFV1_0(w[10], w[55], w[29], pars->GC_11, amp[194]); 
  FFV1_0(w[10], w[78], w[4], pars->GC_11, amp[195]); 
  FFV1_0(w[10], w[54], w[77], pars->GC_11, amp[196]); 
  FFV1_0(w[48], w[55], w[8], pars->GC_2, amp[197]); 
  FFV1_0(w[48], w[65], w[4], pars->GC_11, amp[198]); 
  FFV1_0(w[3], w[65], w[77], pars->GC_11, amp[199]); 
  FFV1_0(w[18], w[55], w[29], pars->GC_11, amp[200]); 
  FFV1_0(w[18], w[78], w[4], pars->GC_11, amp[201]); 
  FFV1_0(w[18], w[54], w[77], pars->GC_11, amp[202]); 
  FFV2_5_0(w[48], w[55], w[17], pars->GC_51, pars->GC_58, amp[203]); 
  FFV1_0(w[3], w[61], w[79], pars->GC_11, amp[204]); 
  FFV1_0(w[3], w[61], w[80], pars->GC_11, amp[205]); 
  FFV1_0(w[3], w[61], w[81], pars->GC_11, amp[206]); 
  FFV1_0(w[10], w[54], w[79], pars->GC_11, amp[207]); 
  FFV1_0(w[10], w[54], w[80], pars->GC_11, amp[208]); 
  FFV1_0(w[10], w[54], w[81], pars->GC_11, amp[209]); 
  FFV1_0(w[3], w[65], w[79], pars->GC_11, amp[210]); 
  FFV1_0(w[3], w[65], w[80], pars->GC_11, amp[211]); 
  FFV1_0(w[3], w[65], w[81], pars->GC_11, amp[212]); 
  FFV1_0(w[18], w[54], w[79], pars->GC_11, amp[213]); 
  FFV1_0(w[18], w[54], w[80], pars->GC_11, amp[214]); 
  FFV1_0(w[18], w[54], w[81], pars->GC_11, amp[215]); 
  FFV1_0(w[15], w[83], w[7], pars->GC_11, amp[216]); 
  FFV1_0(w[13], w[83], w[5], pars->GC_11, amp[217]); 
  FFV1_0(w[57], w[84], w[7], pars->GC_11, amp[218]); 
  FFV1_0(w[13], w[84], w[4], pars->GC_11, amp[219]); 
  FFV1_0(w[57], w[85], w[5], pars->GC_11, amp[220]); 
  FFV1_0(w[15], w[85], w[4], pars->GC_11, amp[221]); 
  FFV1_0(w[21], w[83], w[7], pars->GC_11, amp[222]); 
  FFV1_0(w[20], w[83], w[5], pars->GC_11, amp[223]); 
  FFV1_0(w[59], w[84], w[7], pars->GC_11, amp[224]); 
  FFV1_0(w[20], w[84], w[4], pars->GC_11, amp[225]); 
  FFV1_0(w[59], w[85], w[5], pars->GC_11, amp[226]); 
  FFV1_0(w[21], w[85], w[4], pars->GC_11, amp[227]); 
  FFV1_0(w[87], w[86], w[7], pars->GC_11, amp[228]); 
  FFV1_0(w[63], w[86], w[5], pars->GC_11, amp[229]); 
  FFV1_0(w[64], w[84], w[7], pars->GC_11, amp[230]); 
  FFV1_0(w[64], w[85], w[5], pars->GC_11, amp[231]); 
  FFV1_0(w[63], w[84], w[8], pars->GC_2, amp[232]); 
  FFV1_0(w[87], w[85], w[8], pars->GC_2, amp[233]); 
  FFV1_0(w[87], w[88], w[7], pars->GC_11, amp[234]); 
  FFV1_0(w[63], w[88], w[5], pars->GC_11, amp[235]); 
  FFV1_0(w[66], w[84], w[7], pars->GC_11, amp[236]); 
  FFV1_0(w[66], w[85], w[5], pars->GC_11, amp[237]); 
  FFV2_5_0(w[63], w[84], w[17], pars->GC_51, pars->GC_58, amp[238]); 
  FFV2_5_0(w[87], w[85], w[17], pars->GC_51, pars->GC_58, amp[239]); 
  FFV1_0(w[60], w[86], w[35], pars->GC_11, amp[240]); 
  FFV1_0(w[64], w[82], w[35], pars->GC_11, amp[241]); 
  FFV1_0(w[60], w[88], w[35], pars->GC_11, amp[242]); 
  FFV1_0(w[66], w[82], w[35], pars->GC_11, amp[243]); 
  FFV1_0(w[89], w[86], w[7], pars->GC_11, amp[244]); 
  FFV1_0(w[25], w[86], w[4], pars->GC_11, amp[245]); 
  FFV1_0(w[26], w[83], w[7], pars->GC_11, amp[246]); 
  FFV1_0(w[26], w[85], w[4], pars->GC_11, amp[247]); 
  FFV1_0(w[25], w[83], w[8], pars->GC_2, amp[248]); 
  FFV1_0(w[89], w[85], w[8], pars->GC_2, amp[249]); 
  FFV1_0(w[89], w[88], w[7], pars->GC_11, amp[250]); 
  FFV1_0(w[25], w[88], w[4], pars->GC_11, amp[251]); 
  FFV1_0(w[28], w[83], w[7], pars->GC_11, amp[252]); 
  FFV1_0(w[28], w[85], w[4], pars->GC_11, amp[253]); 
  FFV2_5_0(w[25], w[83], w[17], pars->GC_51, pars->GC_58, amp[254]); 
  FFV2_5_0(w[89], w[85], w[17], pars->GC_51, pars->GC_58, amp[255]); 
  FFV1_0(w[22], w[86], w[68], pars->GC_11, amp[256]); 
  FFV1_0(w[26], w[82], w[68], pars->GC_11, amp[257]); 
  FFV1_0(w[22], w[88], w[68], pars->GC_11, amp[258]); 
  FFV1_0(w[28], w[82], w[68], pars->GC_11, amp[259]); 
  FFV1_0(w[69], w[86], w[5], pars->GC_11, amp[260]); 
  FFV1_0(w[37], w[86], w[4], pars->GC_11, amp[261]); 
  FFV1_0(w[39], w[83], w[5], pars->GC_11, amp[262]); 
  FFV1_0(w[39], w[84], w[4], pars->GC_11, amp[263]); 
  FFV1_0(w[37], w[83], w[8], pars->GC_2, amp[264]); 
  FFV1_0(w[69], w[84], w[8], pars->GC_2, amp[265]); 
  FFV1_0(w[69], w[88], w[5], pars->GC_11, amp[266]); 
  FFV1_0(w[37], w[88], w[4], pars->GC_11, amp[267]); 
  FFV1_0(w[40], w[83], w[5], pars->GC_11, amp[268]); 
  FFV1_0(w[40], w[84], w[4], pars->GC_11, amp[269]); 
  FFV2_5_0(w[37], w[83], w[17], pars->GC_51, pars->GC_58, amp[270]); 
  FFV2_5_0(w[69], w[84], w[17], pars->GC_51, pars->GC_58, amp[271]); 
  FFV1_0(w[36], w[86], w[90], pars->GC_11, amp[272]); 
  FFV1_0(w[39], w[82], w[90], pars->GC_11, amp[273]); 
  FFV1_0(w[36], w[88], w[90], pars->GC_11, amp[274]); 
  FFV1_0(w[40], w[82], w[90], pars->GC_11, amp[275]); 
  FFV1_0(w[91], w[86], w[7], pars->GC_11, amp[276]); 
  FFV1_0(w[3], w[86], w[92], pars->GC_11, amp[277]); 
  FFV1_0(w[10], w[93], w[7], pars->GC_11, amp[278]); 
  FFV1_0(w[10], w[85], w[90], pars->GC_11, amp[279]); 
  FFV1_0(w[10], w[82], w[92], pars->GC_11, amp[280]); 
  FFV1_0(w[91], w[85], w[8], pars->GC_2, amp[281]); 
  FFV1_0(w[91], w[88], w[7], pars->GC_11, amp[282]); 
  FFV1_0(w[3], w[88], w[92], pars->GC_11, amp[283]); 
  FFV1_0(w[18], w[93], w[7], pars->GC_11, amp[284]); 
  FFV1_0(w[18], w[85], w[90], pars->GC_11, amp[285]); 
  FFV1_0(w[18], w[82], w[92], pars->GC_11, amp[286]); 
  FFV2_5_0(w[91], w[85], w[17], pars->GC_51, pars->GC_58, amp[287]); 
  FFV1_0(w[74], w[86], w[5], pars->GC_11, amp[288]); 
  FFV1_0(w[3], w[86], w[94], pars->GC_11, amp[289]); 
  FFV1_0(w[10], w[95], w[5], pars->GC_11, amp[290]); 
  FFV1_0(w[10], w[84], w[68], pars->GC_11, amp[291]); 
  FFV1_0(w[10], w[82], w[94], pars->GC_11, amp[292]); 
  FFV1_0(w[74], w[84], w[8], pars->GC_2, amp[293]); 
  FFV1_0(w[74], w[88], w[5], pars->GC_11, amp[294]); 
  FFV1_0(w[3], w[88], w[94], pars->GC_11, amp[295]); 
  FFV1_0(w[18], w[95], w[5], pars->GC_11, amp[296]); 
  FFV1_0(w[18], w[84], w[68], pars->GC_11, amp[297]); 
  FFV1_0(w[18], w[82], w[94], pars->GC_11, amp[298]); 
  FFV2_5_0(w[74], w[84], w[17], pars->GC_51, pars->GC_58, amp[299]); 
  FFV1_0(w[45], w[86], w[4], pars->GC_11, amp[300]); 
  FFV1_0(w[3], w[86], w[96], pars->GC_11, amp[301]); 
  FFV1_0(w[10], w[83], w[35], pars->GC_11, amp[302]); 
  FFV1_0(w[10], w[97], w[4], pars->GC_11, amp[303]); 
  FFV1_0(w[10], w[82], w[96], pars->GC_11, amp[304]); 
  FFV1_0(w[45], w[83], w[8], pars->GC_2, amp[305]); 
  FFV1_0(w[45], w[88], w[4], pars->GC_11, amp[306]); 
  FFV1_0(w[3], w[88], w[96], pars->GC_11, amp[307]); 
  FFV1_0(w[18], w[83], w[35], pars->GC_11, amp[308]); 
  FFV1_0(w[18], w[97], w[4], pars->GC_11, amp[309]); 
  FFV1_0(w[18], w[82], w[96], pars->GC_11, amp[310]); 
  FFV2_5_0(w[45], w[83], w[17], pars->GC_51, pars->GC_58, amp[311]); 
  FFV1_0(w[3], w[86], w[98], pars->GC_11, amp[312]); 
  FFV1_0(w[3], w[86], w[99], pars->GC_11, amp[313]); 
  FFV1_0(w[3], w[86], w[100], pars->GC_11, amp[314]); 
  FFV1_0(w[10], w[82], w[98], pars->GC_11, amp[315]); 
  FFV1_0(w[10], w[82], w[99], pars->GC_11, amp[316]); 
  FFV1_0(w[10], w[82], w[100], pars->GC_11, amp[317]); 
  FFV1_0(w[3], w[88], w[98], pars->GC_11, amp[318]); 
  FFV1_0(w[3], w[88], w[99], pars->GC_11, amp[319]); 
  FFV1_0(w[3], w[88], w[100], pars->GC_11, amp[320]); 
  FFV1_0(w[18], w[82], w[98], pars->GC_11, amp[321]); 
  FFV1_0(w[18], w[82], w[99], pars->GC_11, amp[322]); 
  FFV1_0(w[18], w[82], w[100], pars->GC_11, amp[323]); 
  FFV1_0(w[15], w[102], w[6], pars->GC_11, amp[324]); 
  FFV1_0(w[12], w[102], w[5], pars->GC_11, amp[325]); 
  FFV1_0(w[57], w[103], w[6], pars->GC_11, amp[326]); 
  FFV1_0(w[12], w[103], w[4], pars->GC_11, amp[327]); 
  FFV1_0(w[57], w[104], w[5], pars->GC_11, amp[328]); 
  FFV1_0(w[15], w[104], w[4], pars->GC_11, amp[329]); 
  FFV1_0(w[21], w[102], w[6], pars->GC_11, amp[330]); 
  FFV1_0(w[19], w[102], w[5], pars->GC_11, amp[331]); 
  FFV1_0(w[59], w[103], w[6], pars->GC_11, amp[332]); 
  FFV1_0(w[19], w[103], w[4], pars->GC_11, amp[333]); 
  FFV1_0(w[59], w[104], w[5], pars->GC_11, amp[334]); 
  FFV1_0(w[21], w[104], w[4], pars->GC_11, amp[335]); 
  FFV1_0(w[87], w[105], w[6], pars->GC_11, amp[336]); 
  FFV1_0(w[62], w[105], w[5], pars->GC_11, amp[337]); 
  FFV1_0(w[64], w[103], w[6], pars->GC_11, amp[338]); 
  FFV1_0(w[64], w[104], w[5], pars->GC_11, amp[339]); 
  FFV1_0(w[62], w[103], w[8], pars->GC_2, amp[340]); 
  FFV1_0(w[87], w[104], w[8], pars->GC_2, amp[341]); 
  FFV1_0(w[87], w[106], w[6], pars->GC_11, amp[342]); 
  FFV1_0(w[62], w[106], w[5], pars->GC_11, amp[343]); 
  FFV1_0(w[66], w[103], w[6], pars->GC_11, amp[344]); 
  FFV1_0(w[66], w[104], w[5], pars->GC_11, amp[345]); 
  FFV2_5_0(w[62], w[103], w[17], pars->GC_51, pars->GC_58, amp[346]); 
  FFV2_5_0(w[87], w[104], w[17], pars->GC_51, pars->GC_58, amp[347]); 
  FFV1_0(w[60], w[105], w[41], pars->GC_11, amp[348]); 
  FFV1_0(w[64], w[101], w[41], pars->GC_11, amp[349]); 
  FFV1_0(w[60], w[106], w[41], pars->GC_11, amp[350]); 
  FFV1_0(w[66], w[101], w[41], pars->GC_11, amp[351]); 
  FFV1_0(w[89], w[105], w[6], pars->GC_11, amp[352]); 
  FFV1_0(w[24], w[105], w[4], pars->GC_11, amp[353]); 
  FFV1_0(w[26], w[102], w[6], pars->GC_11, amp[354]); 
  FFV1_0(w[26], w[104], w[4], pars->GC_11, amp[355]); 
  FFV1_0(w[24], w[102], w[8], pars->GC_2, amp[356]); 
  FFV1_0(w[89], w[104], w[8], pars->GC_2, amp[357]); 
  FFV1_0(w[89], w[106], w[6], pars->GC_11, amp[358]); 
  FFV1_0(w[24], w[106], w[4], pars->GC_11, amp[359]); 
  FFV1_0(w[28], w[102], w[6], pars->GC_11, amp[360]); 
  FFV1_0(w[28], w[104], w[4], pars->GC_11, amp[361]); 
  FFV2_5_0(w[24], w[102], w[17], pars->GC_51, pars->GC_58, amp[362]); 
  FFV2_5_0(w[89], w[104], w[17], pars->GC_51, pars->GC_58, amp[363]); 
  FFV1_0(w[22], w[105], w[70], pars->GC_11, amp[364]); 
  FFV1_0(w[26], w[101], w[70], pars->GC_11, amp[365]); 
  FFV1_0(w[22], w[106], w[70], pars->GC_11, amp[366]); 
  FFV1_0(w[28], w[101], w[70], pars->GC_11, amp[367]); 
  FFV1_0(w[67], w[105], w[5], pars->GC_11, amp[368]); 
  FFV1_0(w[31], w[105], w[4], pars->GC_11, amp[369]); 
  FFV1_0(w[33], w[102], w[5], pars->GC_11, amp[370]); 
  FFV1_0(w[33], w[103], w[4], pars->GC_11, amp[371]); 
  FFV1_0(w[31], w[102], w[8], pars->GC_2, amp[372]); 
  FFV1_0(w[67], w[103], w[8], pars->GC_2, amp[373]); 
  FFV1_0(w[67], w[106], w[5], pars->GC_11, amp[374]); 
  FFV1_0(w[31], w[106], w[4], pars->GC_11, amp[375]); 
  FFV1_0(w[34], w[102], w[5], pars->GC_11, amp[376]); 
  FFV1_0(w[34], w[103], w[4], pars->GC_11, amp[377]); 
  FFV2_5_0(w[31], w[102], w[17], pars->GC_51, pars->GC_58, amp[378]); 
  FFV2_5_0(w[67], w[103], w[17], pars->GC_51, pars->GC_58, amp[379]); 
  FFV1_0(w[30], w[105], w[90], pars->GC_11, amp[380]); 
  FFV1_0(w[33], w[101], w[90], pars->GC_11, amp[381]); 
  FFV1_0(w[30], w[106], w[90], pars->GC_11, amp[382]); 
  FFV1_0(w[34], w[101], w[90], pars->GC_11, amp[383]); 
  FFV1_0(w[91], w[105], w[6], pars->GC_11, amp[384]); 
  FFV1_0(w[3], w[105], w[107], pars->GC_11, amp[385]); 
  FFV1_0(w[10], w[108], w[6], pars->GC_11, amp[386]); 
  FFV1_0(w[10], w[104], w[90], pars->GC_11, amp[387]); 
  FFV1_0(w[10], w[101], w[107], pars->GC_11, amp[388]); 
  FFV1_0(w[91], w[104], w[8], pars->GC_2, amp[389]); 
  FFV1_0(w[91], w[106], w[6], pars->GC_11, amp[390]); 
  FFV1_0(w[3], w[106], w[107], pars->GC_11, amp[391]); 
  FFV1_0(w[18], w[108], w[6], pars->GC_11, amp[392]); 
  FFV1_0(w[18], w[104], w[90], pars->GC_11, amp[393]); 
  FFV1_0(w[18], w[101], w[107], pars->GC_11, amp[394]); 
  FFV2_5_0(w[91], w[104], w[17], pars->GC_51, pars->GC_58, amp[395]); 
  FFV1_0(w[71], w[105], w[5], pars->GC_11, amp[396]); 
  FFV1_0(w[3], w[105], w[109], pars->GC_11, amp[397]); 
  FFV1_0(w[10], w[110], w[5], pars->GC_11, amp[398]); 
  FFV1_0(w[10], w[103], w[70], pars->GC_11, amp[399]); 
  FFV1_0(w[10], w[101], w[109], pars->GC_11, amp[400]); 
  FFV1_0(w[71], w[103], w[8], pars->GC_2, amp[401]); 
  FFV1_0(w[71], w[106], w[5], pars->GC_11, amp[402]); 
  FFV1_0(w[3], w[106], w[109], pars->GC_11, amp[403]); 
  FFV1_0(w[18], w[110], w[5], pars->GC_11, amp[404]); 
  FFV1_0(w[18], w[103], w[70], pars->GC_11, amp[405]); 
  FFV1_0(w[18], w[101], w[109], pars->GC_11, amp[406]); 
  FFV2_5_0(w[71], w[103], w[17], pars->GC_51, pars->GC_58, amp[407]); 
  FFV1_0(w[42], w[105], w[4], pars->GC_11, amp[408]); 
  FFV1_0(w[3], w[105], w[111], pars->GC_11, amp[409]); 
  FFV1_0(w[10], w[102], w[41], pars->GC_11, amp[410]); 
  FFV1_0(w[10], w[112], w[4], pars->GC_11, amp[411]); 
  FFV1_0(w[10], w[101], w[111], pars->GC_11, amp[412]); 
  FFV1_0(w[42], w[102], w[8], pars->GC_2, amp[413]); 
  FFV1_0(w[42], w[106], w[4], pars->GC_11, amp[414]); 
  FFV1_0(w[3], w[106], w[111], pars->GC_11, amp[415]); 
  FFV1_0(w[18], w[102], w[41], pars->GC_11, amp[416]); 
  FFV1_0(w[18], w[112], w[4], pars->GC_11, amp[417]); 
  FFV1_0(w[18], w[101], w[111], pars->GC_11, amp[418]); 
  FFV2_5_0(w[42], w[102], w[17], pars->GC_51, pars->GC_58, amp[419]); 
  FFV1_0(w[3], w[105], w[113], pars->GC_11, amp[420]); 
  FFV1_0(w[3], w[105], w[114], pars->GC_11, amp[421]); 
  FFV1_0(w[3], w[105], w[115], pars->GC_11, amp[422]); 
  FFV1_0(w[10], w[101], w[113], pars->GC_11, amp[423]); 
  FFV1_0(w[10], w[101], w[114], pars->GC_11, amp[424]); 
  FFV1_0(w[10], w[101], w[115], pars->GC_11, amp[425]); 
  FFV1_0(w[3], w[106], w[113], pars->GC_11, amp[426]); 
  FFV1_0(w[3], w[106], w[114], pars->GC_11, amp[427]); 
  FFV1_0(w[3], w[106], w[115], pars->GC_11, amp[428]); 
  FFV1_0(w[18], w[101], w[113], pars->GC_11, amp[429]); 
  FFV1_0(w[18], w[101], w[114], pars->GC_11, amp[430]); 
  FFV1_0(w[18], w[101], w[115], pars->GC_11, amp[431]); 
  FFV1_0(w[87], w[117], w[7], pars->GC_11, amp[432]); 
  FFV1_0(w[87], w[118], w[6], pars->GC_11, amp[433]); 
  FFV1_0(w[62], w[119], w[7], pars->GC_11, amp[434]); 
  FFV1_0(w[62], w[118], w[5], pars->GC_11, amp[435]); 
  FFV1_0(w[63], w[119], w[6], pars->GC_11, amp[436]); 
  FFV1_0(w[63], w[117], w[5], pars->GC_11, amp[437]); 
  FFV1_0(w[87], w[121], w[7], pars->GC_11, amp[438]); 
  FFV1_0(w[87], w[122], w[6], pars->GC_11, amp[439]); 
  FFV1_0(w[62], w[123], w[7], pars->GC_11, amp[440]); 
  FFV1_0(w[62], w[122], w[5], pars->GC_11, amp[441]); 
  FFV1_0(w[63], w[123], w[6], pars->GC_11, amp[442]); 
  FFV1_0(w[63], w[121], w[5], pars->GC_11, amp[443]); 
  FFV1_0(w[124], w[116], w[7], pars->GC_11, amp[444]); 
  FFV1_0(w[63], w[116], w[41], pars->GC_11, amp[445]); 
  FFV1_0(w[60], w[116], w[43], pars->GC_11, amp[446]); 
  FFV1_0(w[64], w[125], w[7], pars->GC_11, amp[447]); 
  FFV1_0(w[64], w[2], w[43], pars->GC_11, amp[448]); 
  FFV1_0(w[63], w[125], w[8], pars->GC_2, amp[449]); 
  FFV1_0(w[124], w[120], w[7], pars->GC_11, amp[450]); 
  FFV1_0(w[63], w[120], w[41], pars->GC_11, amp[451]); 
  FFV1_0(w[60], w[120], w[43], pars->GC_11, amp[452]); 
  FFV1_0(w[66], w[125], w[7], pars->GC_11, amp[453]); 
  FFV1_0(w[66], w[2], w[43], pars->GC_11, amp[454]); 
  FFV2_5_0(w[63], w[125], w[17], pars->GC_51, pars->GC_58, amp[455]); 
  FFV1_0(w[126], w[116], w[6], pars->GC_11, amp[456]); 
  FFV1_0(w[62], w[116], w[35], pars->GC_11, amp[457]); 
  FFV1_0(w[60], w[116], w[46], pars->GC_11, amp[458]); 
  FFV1_0(w[64], w[127], w[6], pars->GC_11, amp[459]); 
  FFV1_0(w[64], w[2], w[46], pars->GC_11, amp[460]); 
  FFV1_0(w[62], w[127], w[8], pars->GC_2, amp[461]); 
  FFV1_0(w[126], w[120], w[6], pars->GC_11, amp[462]); 
  FFV1_0(w[62], w[120], w[35], pars->GC_11, amp[463]); 
  FFV1_0(w[60], w[120], w[46], pars->GC_11, amp[464]); 
  FFV1_0(w[66], w[127], w[6], pars->GC_11, amp[465]); 
  FFV1_0(w[66], w[2], w[46], pars->GC_11, amp[466]); 
  FFV2_5_0(w[62], w[127], w[17], pars->GC_51, pars->GC_58, amp[467]); 
  FFV1_0(w[87], w[116], w[29], pars->GC_11, amp[468]); 
  FFV1_0(w[128], w[116], w[5], pars->GC_11, amp[469]); 
  FFV1_0(w[60], w[116], w[49], pars->GC_11, amp[470]); 
  FFV1_0(w[64], w[129], w[5], pars->GC_11, amp[471]); 
  FFV1_0(w[64], w[2], w[49], pars->GC_11, amp[472]); 
  FFV1_0(w[87], w[129], w[8], pars->GC_2, amp[473]); 
  FFV1_0(w[87], w[120], w[29], pars->GC_11, amp[474]); 
  FFV1_0(w[128], w[120], w[5], pars->GC_11, amp[475]); 
  FFV1_0(w[60], w[120], w[49], pars->GC_11, amp[476]); 
  FFV1_0(w[66], w[129], w[5], pars->GC_11, amp[477]); 
  FFV1_0(w[66], w[2], w[49], pars->GC_11, amp[478]); 
  FFV2_5_0(w[87], w[129], w[17], pars->GC_51, pars->GC_58, amp[479]); 
  FFV1_0(w[60], w[116], w[51], pars->GC_11, amp[480]); 
  FFV1_0(w[60], w[116], w[52], pars->GC_11, amp[481]); 
  FFV1_0(w[60], w[116], w[53], pars->GC_11, amp[482]); 
  FFV1_0(w[64], w[2], w[51], pars->GC_11, amp[483]); 
  FFV1_0(w[64], w[2], w[52], pars->GC_11, amp[484]); 
  FFV1_0(w[64], w[2], w[53], pars->GC_11, amp[485]); 
  FFV1_0(w[60], w[120], w[51], pars->GC_11, amp[486]); 
  FFV1_0(w[60], w[120], w[52], pars->GC_11, amp[487]); 
  FFV1_0(w[60], w[120], w[53], pars->GC_11, amp[488]); 
  FFV1_0(w[66], w[2], w[51], pars->GC_11, amp[489]); 
  FFV1_0(w[66], w[2], w[52], pars->GC_11, amp[490]); 
  FFV1_0(w[66], w[2], w[53], pars->GC_11, amp[491]); 
  FFV1_0(w[89], w[117], w[7], pars->GC_11, amp[492]); 
  FFV1_0(w[89], w[118], w[6], pars->GC_11, amp[493]); 
  FFV1_0(w[24], w[130], w[7], pars->GC_11, amp[494]); 
  FFV1_0(w[24], w[118], w[4], pars->GC_11, amp[495]); 
  FFV1_0(w[25], w[130], w[6], pars->GC_11, amp[496]); 
  FFV1_0(w[25], w[117], w[4], pars->GC_11, amp[497]); 
  FFV1_0(w[89], w[121], w[7], pars->GC_11, amp[498]); 
  FFV1_0(w[89], w[122], w[6], pars->GC_11, amp[499]); 
  FFV1_0(w[24], w[131], w[7], pars->GC_11, amp[500]); 
  FFV1_0(w[24], w[122], w[4], pars->GC_11, amp[501]); 
  FFV1_0(w[25], w[131], w[6], pars->GC_11, amp[502]); 
  FFV1_0(w[25], w[121], w[4], pars->GC_11, amp[503]); 
  FFV1_0(w[132], w[116], w[7], pars->GC_11, amp[504]); 
  FFV1_0(w[25], w[116], w[70], pars->GC_11, amp[505]); 
  FFV1_0(w[22], w[116], w[72], pars->GC_11, amp[506]); 
  FFV1_0(w[26], w[133], w[7], pars->GC_11, amp[507]); 
  FFV1_0(w[26], w[2], w[72], pars->GC_11, amp[508]); 
  FFV1_0(w[25], w[133], w[8], pars->GC_2, amp[509]); 
  FFV1_0(w[132], w[120], w[7], pars->GC_11, amp[510]); 
  FFV1_0(w[25], w[120], w[70], pars->GC_11, amp[511]); 
  FFV1_0(w[22], w[120], w[72], pars->GC_11, amp[512]); 
  FFV1_0(w[28], w[133], w[7], pars->GC_11, amp[513]); 
  FFV1_0(w[28], w[2], w[72], pars->GC_11, amp[514]); 
  FFV2_5_0(w[25], w[133], w[17], pars->GC_51, pars->GC_58, amp[515]); 
  FFV1_0(w[134], w[116], w[6], pars->GC_11, amp[516]); 
  FFV1_0(w[24], w[116], w[68], pars->GC_11, amp[517]); 
  FFV1_0(w[22], w[116], w[75], pars->GC_11, amp[518]); 
  FFV1_0(w[26], w[135], w[6], pars->GC_11, amp[519]); 
  FFV1_0(w[26], w[2], w[75], pars->GC_11, amp[520]); 
  FFV1_0(w[24], w[135], w[8], pars->GC_2, amp[521]); 
  FFV1_0(w[134], w[120], w[6], pars->GC_11, amp[522]); 
  FFV1_0(w[24], w[120], w[68], pars->GC_11, amp[523]); 
  FFV1_0(w[22], w[120], w[75], pars->GC_11, amp[524]); 
  FFV1_0(w[28], w[135], w[6], pars->GC_11, amp[525]); 
  FFV1_0(w[28], w[2], w[75], pars->GC_11, amp[526]); 
  FFV2_5_0(w[24], w[135], w[17], pars->GC_51, pars->GC_58, amp[527]); 
  FFV1_0(w[89], w[116], w[29], pars->GC_11, amp[528]); 
  FFV1_0(w[136], w[116], w[4], pars->GC_11, amp[529]); 
  FFV1_0(w[22], w[116], w[77], pars->GC_11, amp[530]); 
  FFV1_0(w[26], w[129], w[4], pars->GC_11, amp[531]); 
  FFV1_0(w[26], w[2], w[77], pars->GC_11, amp[532]); 
  FFV1_0(w[89], w[129], w[8], pars->GC_2, amp[533]); 
  FFV1_0(w[89], w[120], w[29], pars->GC_11, amp[534]); 
  FFV1_0(w[136], w[120], w[4], pars->GC_11, amp[535]); 
  FFV1_0(w[22], w[120], w[77], pars->GC_11, amp[536]); 
  FFV1_0(w[28], w[129], w[4], pars->GC_11, amp[537]); 
  FFV1_0(w[28], w[2], w[77], pars->GC_11, amp[538]); 
  FFV2_5_0(w[89], w[129], w[17], pars->GC_51, pars->GC_58, amp[539]); 
  FFV1_0(w[22], w[116], w[79], pars->GC_11, amp[540]); 
  FFV1_0(w[22], w[116], w[80], pars->GC_11, amp[541]); 
  FFV1_0(w[22], w[116], w[81], pars->GC_11, amp[542]); 
  FFV1_0(w[26], w[2], w[79], pars->GC_11, amp[543]); 
  FFV1_0(w[26], w[2], w[80], pars->GC_11, amp[544]); 
  FFV1_0(w[26], w[2], w[81], pars->GC_11, amp[545]); 
  FFV1_0(w[22], w[120], w[79], pars->GC_11, amp[546]); 
  FFV1_0(w[22], w[120], w[80], pars->GC_11, amp[547]); 
  FFV1_0(w[22], w[120], w[81], pars->GC_11, amp[548]); 
  FFV1_0(w[28], w[2], w[79], pars->GC_11, amp[549]); 
  FFV1_0(w[28], w[2], w[80], pars->GC_11, amp[550]); 
  FFV1_0(w[28], w[2], w[81], pars->GC_11, amp[551]); 
  FFV1_0(w[67], w[119], w[7], pars->GC_11, amp[552]); 
  FFV1_0(w[67], w[118], w[5], pars->GC_11, amp[553]); 
  FFV1_0(w[31], w[130], w[7], pars->GC_11, amp[554]); 
  FFV1_0(w[31], w[118], w[4], pars->GC_11, amp[555]); 
  FFV1_0(w[32], w[130], w[5], pars->GC_11, amp[556]); 
  FFV1_0(w[32], w[119], w[4], pars->GC_11, amp[557]); 
  FFV1_0(w[67], w[123], w[7], pars->GC_11, amp[558]); 
  FFV1_0(w[67], w[122], w[5], pars->GC_11, amp[559]); 
  FFV1_0(w[31], w[131], w[7], pars->GC_11, amp[560]); 
  FFV1_0(w[31], w[122], w[4], pars->GC_11, amp[561]); 
  FFV1_0(w[32], w[131], w[5], pars->GC_11, amp[562]); 
  FFV1_0(w[32], w[123], w[4], pars->GC_11, amp[563]); 
  FFV1_0(w[137], w[116], w[7], pars->GC_11, amp[564]); 
  FFV1_0(w[32], w[116], w[90], pars->GC_11, amp[565]); 
  FFV1_0(w[30], w[116], w[92], pars->GC_11, amp[566]); 
  FFV1_0(w[33], w[138], w[7], pars->GC_11, amp[567]); 
  FFV1_0(w[33], w[2], w[92], pars->GC_11, amp[568]); 
  FFV1_0(w[32], w[138], w[8], pars->GC_2, amp[569]); 
  FFV1_0(w[137], w[120], w[7], pars->GC_11, amp[570]); 
  FFV1_0(w[32], w[120], w[90], pars->GC_11, amp[571]); 
  FFV1_0(w[30], w[120], w[92], pars->GC_11, amp[572]); 
  FFV1_0(w[34], w[138], w[7], pars->GC_11, amp[573]); 
  FFV1_0(w[34], w[2], w[92], pars->GC_11, amp[574]); 
  FFV2_5_0(w[32], w[138], w[17], pars->GC_51, pars->GC_58, amp[575]); 
  FFV1_0(w[139], w[116], w[5], pars->GC_11, amp[576]); 
  FFV1_0(w[31], w[116], w[68], pars->GC_11, amp[577]); 
  FFV1_0(w[30], w[116], w[94], pars->GC_11, amp[578]); 
  FFV1_0(w[33], w[135], w[5], pars->GC_11, amp[579]); 
  FFV1_0(w[33], w[2], w[94], pars->GC_11, amp[580]); 
  FFV1_0(w[31], w[135], w[8], pars->GC_2, amp[581]); 
  FFV1_0(w[139], w[120], w[5], pars->GC_11, amp[582]); 
  FFV1_0(w[31], w[120], w[68], pars->GC_11, amp[583]); 
  FFV1_0(w[30], w[120], w[94], pars->GC_11, amp[584]); 
  FFV1_0(w[34], w[135], w[5], pars->GC_11, amp[585]); 
  FFV1_0(w[34], w[2], w[94], pars->GC_11, amp[586]); 
  FFV2_5_0(w[31], w[135], w[17], pars->GC_51, pars->GC_58, amp[587]); 
  FFV1_0(w[67], w[116], w[35], pars->GC_11, amp[588]); 
  FFV1_0(w[140], w[116], w[4], pars->GC_11, amp[589]); 
  FFV1_0(w[30], w[116], w[96], pars->GC_11, amp[590]); 
  FFV1_0(w[33], w[127], w[4], pars->GC_11, amp[591]); 
  FFV1_0(w[33], w[2], w[96], pars->GC_11, amp[592]); 
  FFV1_0(w[67], w[127], w[8], pars->GC_2, amp[593]); 
  FFV1_0(w[67], w[120], w[35], pars->GC_11, amp[594]); 
  FFV1_0(w[140], w[120], w[4], pars->GC_11, amp[595]); 
  FFV1_0(w[30], w[120], w[96], pars->GC_11, amp[596]); 
  FFV1_0(w[34], w[127], w[4], pars->GC_11, amp[597]); 
  FFV1_0(w[34], w[2], w[96], pars->GC_11, amp[598]); 
  FFV2_5_0(w[67], w[127], w[17], pars->GC_51, pars->GC_58, amp[599]); 
  FFV1_0(w[30], w[116], w[98], pars->GC_11, amp[600]); 
  FFV1_0(w[30], w[116], w[99], pars->GC_11, amp[601]); 
  FFV1_0(w[30], w[116], w[100], pars->GC_11, amp[602]); 
  FFV1_0(w[33], w[2], w[98], pars->GC_11, amp[603]); 
  FFV1_0(w[33], w[2], w[99], pars->GC_11, amp[604]); 
  FFV1_0(w[33], w[2], w[100], pars->GC_11, amp[605]); 
  FFV1_0(w[30], w[120], w[98], pars->GC_11, amp[606]); 
  FFV1_0(w[30], w[120], w[99], pars->GC_11, amp[607]); 
  FFV1_0(w[30], w[120], w[100], pars->GC_11, amp[608]); 
  FFV1_0(w[34], w[2], w[98], pars->GC_11, amp[609]); 
  FFV1_0(w[34], w[2], w[99], pars->GC_11, amp[610]); 
  FFV1_0(w[34], w[2], w[100], pars->GC_11, amp[611]); 
  FFV1_0(w[69], w[119], w[6], pars->GC_11, amp[612]); 
  FFV1_0(w[69], w[117], w[5], pars->GC_11, amp[613]); 
  FFV1_0(w[37], w[130], w[6], pars->GC_11, amp[614]); 
  FFV1_0(w[37], w[117], w[4], pars->GC_11, amp[615]); 
  FFV1_0(w[38], w[130], w[5], pars->GC_11, amp[616]); 
  FFV1_0(w[38], w[119], w[4], pars->GC_11, amp[617]); 
  FFV1_0(w[69], w[123], w[6], pars->GC_11, amp[618]); 
  FFV1_0(w[69], w[121], w[5], pars->GC_11, amp[619]); 
  FFV1_0(w[37], w[131], w[6], pars->GC_11, amp[620]); 
  FFV1_0(w[37], w[121], w[4], pars->GC_11, amp[621]); 
  FFV1_0(w[38], w[131], w[5], pars->GC_11, amp[622]); 
  FFV1_0(w[38], w[123], w[4], pars->GC_11, amp[623]); 
  FFV1_0(w[141], w[116], w[6], pars->GC_11, amp[624]); 
  FFV1_0(w[38], w[116], w[90], pars->GC_11, amp[625]); 
  FFV1_0(w[36], w[116], w[107], pars->GC_11, amp[626]); 
  FFV1_0(w[39], w[138], w[6], pars->GC_11, amp[627]); 
  FFV1_0(w[39], w[2], w[107], pars->GC_11, amp[628]); 
  FFV1_0(w[38], w[138], w[8], pars->GC_2, amp[629]); 
  FFV1_0(w[141], w[120], w[6], pars->GC_11, amp[630]); 
  FFV1_0(w[38], w[120], w[90], pars->GC_11, amp[631]); 
  FFV1_0(w[36], w[120], w[107], pars->GC_11, amp[632]); 
  FFV1_0(w[40], w[138], w[6], pars->GC_11, amp[633]); 
  FFV1_0(w[40], w[2], w[107], pars->GC_11, amp[634]); 
  FFV2_5_0(w[38], w[138], w[17], pars->GC_51, pars->GC_58, amp[635]); 
  FFV1_0(w[142], w[116], w[5], pars->GC_11, amp[636]); 
  FFV1_0(w[37], w[116], w[70], pars->GC_11, amp[637]); 
  FFV1_0(w[36], w[116], w[109], pars->GC_11, amp[638]); 
  FFV1_0(w[39], w[133], w[5], pars->GC_11, amp[639]); 
  FFV1_0(w[39], w[2], w[109], pars->GC_11, amp[640]); 
  FFV1_0(w[37], w[133], w[8], pars->GC_2, amp[641]); 
  FFV1_0(w[142], w[120], w[5], pars->GC_11, amp[642]); 
  FFV1_0(w[37], w[120], w[70], pars->GC_11, amp[643]); 
  FFV1_0(w[36], w[120], w[109], pars->GC_11, amp[644]); 
  FFV1_0(w[40], w[133], w[5], pars->GC_11, amp[645]); 
  FFV1_0(w[40], w[2], w[109], pars->GC_11, amp[646]); 
  FFV2_5_0(w[37], w[133], w[17], pars->GC_51, pars->GC_58, amp[647]); 
  FFV1_0(w[69], w[116], w[41], pars->GC_11, amp[648]); 
  FFV1_0(w[143], w[116], w[4], pars->GC_11, amp[649]); 
  FFV1_0(w[36], w[116], w[111], pars->GC_11, amp[650]); 
  FFV1_0(w[39], w[125], w[4], pars->GC_11, amp[651]); 
  FFV1_0(w[39], w[2], w[111], pars->GC_11, amp[652]); 
  FFV1_0(w[69], w[125], w[8], pars->GC_2, amp[653]); 
  FFV1_0(w[69], w[120], w[41], pars->GC_11, amp[654]); 
  FFV1_0(w[143], w[120], w[4], pars->GC_11, amp[655]); 
  FFV1_0(w[36], w[120], w[111], pars->GC_11, amp[656]); 
  FFV1_0(w[40], w[125], w[4], pars->GC_11, amp[657]); 
  FFV1_0(w[40], w[2], w[111], pars->GC_11, amp[658]); 
  FFV2_5_0(w[69], w[125], w[17], pars->GC_51, pars->GC_58, amp[659]); 
  FFV1_0(w[36], w[116], w[113], pars->GC_11, amp[660]); 
  FFV1_0(w[36], w[116], w[114], pars->GC_11, amp[661]); 
  FFV1_0(w[36], w[116], w[115], pars->GC_11, amp[662]); 
  FFV1_0(w[39], w[2], w[113], pars->GC_11, amp[663]); 
  FFV1_0(w[39], w[2], w[114], pars->GC_11, amp[664]); 
  FFV1_0(w[39], w[2], w[115], pars->GC_11, amp[665]); 
  FFV1_0(w[36], w[120], w[113], pars->GC_11, amp[666]); 
  FFV1_0(w[36], w[120], w[114], pars->GC_11, amp[667]); 
  FFV1_0(w[36], w[120], w[115], pars->GC_11, amp[668]); 
  FFV1_0(w[40], w[2], w[113], pars->GC_11, amp[669]); 
  FFV1_0(w[40], w[2], w[114], pars->GC_11, amp[670]); 
  FFV1_0(w[40], w[2], w[115], pars->GC_11, amp[671]); 
  FFV1_0(w[91], w[117], w[7], pars->GC_11, amp[672]); 
  FFV1_0(w[91], w[118], w[6], pars->GC_11, amp[673]); 
  VVV1_0(w[107], w[7], w[144], pars->GC_10, amp[674]); 
  FFV1_0(w[3], w[118], w[107], pars->GC_11, amp[675]); 
  VVV1_0(w[92], w[6], w[144], pars->GC_10, amp[676]); 
  FFV1_0(w[3], w[117], w[92], pars->GC_11, amp[677]); 
  FFV1_0(w[3], w[116], w[145], pars->GC_11, amp[678]); 
  FFV1_0(w[3], w[116], w[146], pars->GC_11, amp[679]); 
  FFV1_0(w[3], w[116], w[147], pars->GC_11, amp[680]); 
  FFV1_0(w[12], w[138], w[7], pars->GC_11, amp[681]); 
  FFV1_0(w[13], w[138], w[6], pars->GC_11, amp[682]); 
  VVV1_0(w[107], w[7], w[148], pars->GC_10, amp[683]); 
  FFV1_0(w[13], w[2], w[107], pars->GC_11, amp[684]); 
  VVV1_0(w[92], w[6], w[148], pars->GC_10, amp[685]); 
  FFV1_0(w[12], w[2], w[92], pars->GC_11, amp[686]); 
  FFV1_0(w[10], w[2], w[145], pars->GC_11, amp[687]); 
  FFV1_0(w[10], w[2], w[146], pars->GC_11, amp[688]); 
  FFV1_0(w[10], w[2], w[147], pars->GC_11, amp[689]); 
  FFV1_0(w[91], w[121], w[7], pars->GC_11, amp[690]); 
  FFV1_0(w[91], w[122], w[6], pars->GC_11, amp[691]); 
  VVV1_0(w[107], w[7], w[149], pars->GC_10, amp[692]); 
  FFV1_0(w[3], w[122], w[107], pars->GC_11, amp[693]); 
  VVV1_0(w[92], w[6], w[149], pars->GC_10, amp[694]); 
  FFV1_0(w[3], w[121], w[92], pars->GC_11, amp[695]); 
  FFV1_0(w[3], w[120], w[145], pars->GC_11, amp[696]); 
  FFV1_0(w[3], w[120], w[146], pars->GC_11, amp[697]); 
  FFV1_0(w[3], w[120], w[147], pars->GC_11, amp[698]); 
  FFV1_0(w[19], w[138], w[7], pars->GC_11, amp[699]); 
  FFV1_0(w[20], w[138], w[6], pars->GC_11, amp[700]); 
  VVV1_0(w[107], w[7], w[150], pars->GC_10, amp[701]); 
  FFV1_0(w[20], w[2], w[107], pars->GC_11, amp[702]); 
  VVV1_0(w[92], w[6], w[150], pars->GC_10, amp[703]); 
  FFV1_0(w[19], w[2], w[92], pars->GC_11, amp[704]); 
  FFV1_0(w[18], w[2], w[145], pars->GC_11, amp[705]); 
  FFV1_0(w[18], w[2], w[146], pars->GC_11, amp[706]); 
  FFV1_0(w[18], w[2], w[147], pars->GC_11, amp[707]); 
  FFV1_0(w[91], w[116], w[29], pars->GC_11, amp[708]); 
  FFV1_0(w[48], w[116], w[90], pars->GC_11, amp[709]); 
  FFV1_0(w[3], w[116], w[151], pars->GC_11, amp[710]); 
  FFV1_0(w[10], w[138], w[29], pars->GC_11, amp[711]); 
  FFV1_0(w[10], w[129], w[90], pars->GC_11, amp[712]); 
  FFV1_0(w[10], w[2], w[151], pars->GC_11, amp[713]); 
  FFV1_0(w[48], w[138], w[8], pars->GC_2, amp[714]); 
  FFV1_0(w[91], w[129], w[8], pars->GC_2, amp[715]); 
  FFV1_0(w[91], w[120], w[29], pars->GC_11, amp[716]); 
  FFV1_0(w[48], w[120], w[90], pars->GC_11, amp[717]); 
  FFV1_0(w[3], w[120], w[151], pars->GC_11, amp[718]); 
  FFV1_0(w[18], w[138], w[29], pars->GC_11, amp[719]); 
  FFV1_0(w[18], w[129], w[90], pars->GC_11, amp[720]); 
  FFV1_0(w[18], w[2], w[151], pars->GC_11, amp[721]); 
  FFV2_5_0(w[48], w[138], w[17], pars->GC_51, pars->GC_58, amp[722]); 
  FFV2_5_0(w[91], w[129], w[17], pars->GC_51, pars->GC_58, amp[723]); 
  FFV1_0(w[71], w[119], w[7], pars->GC_11, amp[724]); 
  FFV1_0(w[71], w[118], w[5], pars->GC_11, amp[725]); 
  VVV1_0(w[109], w[7], w[144], pars->GC_10, amp[726]); 
  FFV1_0(w[3], w[118], w[109], pars->GC_11, amp[727]); 
  VVV1_0(w[72], w[5], w[144], pars->GC_10, amp[728]); 
  FFV1_0(w[3], w[119], w[72], pars->GC_11, amp[729]); 
  FFV1_0(w[3], w[116], w[152], pars->GC_11, amp[730]); 
  FFV1_0(w[3], w[116], w[153], pars->GC_11, amp[731]); 
  FFV1_0(w[3], w[116], w[154], pars->GC_11, amp[732]); 
  FFV1_0(w[15], w[133], w[7], pars->GC_11, amp[733]); 
  FFV1_0(w[13], w[133], w[5], pars->GC_11, amp[734]); 
  VVV1_0(w[109], w[7], w[148], pars->GC_10, amp[735]); 
  FFV1_0(w[13], w[2], w[109], pars->GC_11, amp[736]); 
  VVV1_0(w[72], w[5], w[148], pars->GC_10, amp[737]); 
  FFV1_0(w[15], w[2], w[72], pars->GC_11, amp[738]); 
  FFV1_0(w[10], w[2], w[152], pars->GC_11, amp[739]); 
  FFV1_0(w[10], w[2], w[153], pars->GC_11, amp[740]); 
  FFV1_0(w[10], w[2], w[154], pars->GC_11, amp[741]); 
  FFV1_0(w[71], w[123], w[7], pars->GC_11, amp[742]); 
  FFV1_0(w[71], w[122], w[5], pars->GC_11, amp[743]); 
  VVV1_0(w[109], w[7], w[149], pars->GC_10, amp[744]); 
  FFV1_0(w[3], w[122], w[109], pars->GC_11, amp[745]); 
  VVV1_0(w[72], w[5], w[149], pars->GC_10, amp[746]); 
  FFV1_0(w[3], w[123], w[72], pars->GC_11, amp[747]); 
  FFV1_0(w[3], w[120], w[152], pars->GC_11, amp[748]); 
  FFV1_0(w[3], w[120], w[153], pars->GC_11, amp[749]); 
  FFV1_0(w[3], w[120], w[154], pars->GC_11, amp[750]); 
  FFV1_0(w[21], w[133], w[7], pars->GC_11, amp[751]); 
  FFV1_0(w[20], w[133], w[5], pars->GC_11, amp[752]); 
  VVV1_0(w[109], w[7], w[150], pars->GC_10, amp[753]); 
  FFV1_0(w[20], w[2], w[109], pars->GC_11, amp[754]); 
  VVV1_0(w[72], w[5], w[150], pars->GC_10, amp[755]); 
  FFV1_0(w[21], w[2], w[72], pars->GC_11, amp[756]); 
  FFV1_0(w[18], w[2], w[152], pars->GC_11, amp[757]); 
  FFV1_0(w[18], w[2], w[153], pars->GC_11, amp[758]); 
  FFV1_0(w[18], w[2], w[154], pars->GC_11, amp[759]); 
  FFV1_0(w[71], w[116], w[35], pars->GC_11, amp[760]); 
  FFV1_0(w[45], w[116], w[70], pars->GC_11, amp[761]); 
  FFV1_0(w[3], w[116], w[155], pars->GC_11, amp[762]); 
  FFV1_0(w[10], w[133], w[35], pars->GC_11, amp[763]); 
  FFV1_0(w[10], w[127], w[70], pars->GC_11, amp[764]); 
  FFV1_0(w[10], w[2], w[155], pars->GC_11, amp[765]); 
  FFV1_0(w[45], w[133], w[8], pars->GC_2, amp[766]); 
  FFV1_0(w[71], w[127], w[8], pars->GC_2, amp[767]); 
  FFV1_0(w[71], w[120], w[35], pars->GC_11, amp[768]); 
  FFV1_0(w[45], w[120], w[70], pars->GC_11, amp[769]); 
  FFV1_0(w[3], w[120], w[155], pars->GC_11, amp[770]); 
  FFV1_0(w[18], w[133], w[35], pars->GC_11, amp[771]); 
  FFV1_0(w[18], w[127], w[70], pars->GC_11, amp[772]); 
  FFV1_0(w[18], w[2], w[155], pars->GC_11, amp[773]); 
  FFV2_5_0(w[45], w[133], w[17], pars->GC_51, pars->GC_58, amp[774]); 
  FFV2_5_0(w[71], w[127], w[17], pars->GC_51, pars->GC_58, amp[775]); 
  FFV1_0(w[74], w[119], w[6], pars->GC_11, amp[776]); 
  FFV1_0(w[74], w[117], w[5], pars->GC_11, amp[777]); 
  VVV1_0(w[94], w[6], w[144], pars->GC_10, amp[778]); 
  FFV1_0(w[3], w[117], w[94], pars->GC_11, amp[779]); 
  VVV1_0(w[75], w[5], w[144], pars->GC_10, amp[780]); 
  FFV1_0(w[3], w[119], w[75], pars->GC_11, amp[781]); 
  FFV1_0(w[3], w[116], w[156], pars->GC_11, amp[782]); 
  FFV1_0(w[3], w[116], w[157], pars->GC_11, amp[783]); 
  FFV1_0(w[3], w[116], w[158], pars->GC_11, amp[784]); 
  FFV1_0(w[15], w[135], w[6], pars->GC_11, amp[785]); 
  FFV1_0(w[12], w[135], w[5], pars->GC_11, amp[786]); 
  VVV1_0(w[94], w[6], w[148], pars->GC_10, amp[787]); 
  FFV1_0(w[12], w[2], w[94], pars->GC_11, amp[788]); 
  VVV1_0(w[75], w[5], w[148], pars->GC_10, amp[789]); 
  FFV1_0(w[15], w[2], w[75], pars->GC_11, amp[790]); 
  FFV1_0(w[10], w[2], w[156], pars->GC_11, amp[791]); 
  FFV1_0(w[10], w[2], w[157], pars->GC_11, amp[792]); 
  FFV1_0(w[10], w[2], w[158], pars->GC_11, amp[793]); 
  FFV1_0(w[74], w[123], w[6], pars->GC_11, amp[794]); 
  FFV1_0(w[74], w[121], w[5], pars->GC_11, amp[795]); 
  VVV1_0(w[94], w[6], w[149], pars->GC_10, amp[796]); 
  FFV1_0(w[3], w[121], w[94], pars->GC_11, amp[797]); 
  VVV1_0(w[75], w[5], w[149], pars->GC_10, amp[798]); 
  FFV1_0(w[3], w[123], w[75], pars->GC_11, amp[799]); 
  FFV1_0(w[3], w[120], w[156], pars->GC_11, amp[800]); 
  FFV1_0(w[3], w[120], w[157], pars->GC_11, amp[801]); 
  FFV1_0(w[3], w[120], w[158], pars->GC_11, amp[802]); 
  FFV1_0(w[21], w[135], w[6], pars->GC_11, amp[803]); 
  FFV1_0(w[19], w[135], w[5], pars->GC_11, amp[804]); 
  VVV1_0(w[94], w[6], w[150], pars->GC_10, amp[805]); 
  FFV1_0(w[19], w[2], w[94], pars->GC_11, amp[806]); 
  VVV1_0(w[75], w[5], w[150], pars->GC_10, amp[807]); 
  FFV1_0(w[21], w[2], w[75], pars->GC_11, amp[808]); 
  FFV1_0(w[18], w[2], w[156], pars->GC_11, amp[809]); 
  FFV1_0(w[18], w[2], w[157], pars->GC_11, amp[810]); 
  FFV1_0(w[18], w[2], w[158], pars->GC_11, amp[811]); 
  FFV1_0(w[74], w[116], w[41], pars->GC_11, amp[812]); 
  FFV1_0(w[42], w[116], w[68], pars->GC_11, amp[813]); 
  FFV1_0(w[3], w[116], w[159], pars->GC_11, amp[814]); 
  FFV1_0(w[10], w[135], w[41], pars->GC_11, amp[815]); 
  FFV1_0(w[10], w[125], w[68], pars->GC_11, amp[816]); 
  FFV1_0(w[10], w[2], w[159], pars->GC_11, amp[817]); 
  FFV1_0(w[42], w[135], w[8], pars->GC_2, amp[818]); 
  FFV1_0(w[74], w[125], w[8], pars->GC_2, amp[819]); 
  FFV1_0(w[74], w[120], w[41], pars->GC_11, amp[820]); 
  FFV1_0(w[42], w[120], w[68], pars->GC_11, amp[821]); 
  FFV1_0(w[3], w[120], w[159], pars->GC_11, amp[822]); 
  FFV1_0(w[18], w[135], w[41], pars->GC_11, amp[823]); 
  FFV1_0(w[18], w[125], w[68], pars->GC_11, amp[824]); 
  FFV1_0(w[18], w[2], w[159], pars->GC_11, amp[825]); 
  FFV2_5_0(w[42], w[135], w[17], pars->GC_51, pars->GC_58, amp[826]); 
  FFV2_5_0(w[74], w[125], w[17], pars->GC_51, pars->GC_58, amp[827]); 
  FFV1_0(w[42], w[130], w[7], pars->GC_11, amp[828]); 
  FFV1_0(w[42], w[118], w[4], pars->GC_11, amp[829]); 
  VVV1_0(w[111], w[7], w[144], pars->GC_10, amp[830]); 
  FFV1_0(w[3], w[118], w[111], pars->GC_11, amp[831]); 
  VVV1_0(w[4], w[43], w[144], pars->GC_10, amp[832]); 
  FFV1_0(w[3], w[130], w[43], pars->GC_11, amp[833]); 
  FFV1_0(w[3], w[116], w[160], pars->GC_11, amp[834]); 
  FFV1_0(w[3], w[116], w[161], pars->GC_11, amp[835]); 
  FFV1_0(w[3], w[116], w[162], pars->GC_11, amp[836]); 
  FFV1_0(w[57], w[125], w[7], pars->GC_11, amp[837]); 
  FFV1_0(w[13], w[125], w[4], pars->GC_11, amp[838]); 
  VVV1_0(w[111], w[7], w[148], pars->GC_10, amp[839]); 
  FFV1_0(w[13], w[2], w[111], pars->GC_11, amp[840]); 
  VVV1_0(w[4], w[43], w[148], pars->GC_10, amp[841]); 
  FFV1_0(w[57], w[2], w[43], pars->GC_11, amp[842]); 
  FFV1_0(w[10], w[2], w[160], pars->GC_11, amp[843]); 
  FFV1_0(w[10], w[2], w[161], pars->GC_11, amp[844]); 
  FFV1_0(w[10], w[2], w[162], pars->GC_11, amp[845]); 
  FFV1_0(w[42], w[131], w[7], pars->GC_11, amp[846]); 
  FFV1_0(w[42], w[122], w[4], pars->GC_11, amp[847]); 
  VVV1_0(w[111], w[7], w[149], pars->GC_10, amp[848]); 
  FFV1_0(w[3], w[122], w[111], pars->GC_11, amp[849]); 
  VVV1_0(w[4], w[43], w[149], pars->GC_10, amp[850]); 
  FFV1_0(w[3], w[131], w[43], pars->GC_11, amp[851]); 
  FFV1_0(w[3], w[120], w[160], pars->GC_11, amp[852]); 
  FFV1_0(w[3], w[120], w[161], pars->GC_11, amp[853]); 
  FFV1_0(w[3], w[120], w[162], pars->GC_11, amp[854]); 
  FFV1_0(w[59], w[125], w[7], pars->GC_11, amp[855]); 
  FFV1_0(w[20], w[125], w[4], pars->GC_11, amp[856]); 
  VVV1_0(w[111], w[7], w[150], pars->GC_10, amp[857]); 
  FFV1_0(w[20], w[2], w[111], pars->GC_11, amp[858]); 
  VVV1_0(w[4], w[43], w[150], pars->GC_10, amp[859]); 
  FFV1_0(w[59], w[2], w[43], pars->GC_11, amp[860]); 
  FFV1_0(w[18], w[2], w[160], pars->GC_11, amp[861]); 
  FFV1_0(w[18], w[2], w[161], pars->GC_11, amp[862]); 
  FFV1_0(w[18], w[2], w[162], pars->GC_11, amp[863]); 
  FFV1_0(w[45], w[130], w[6], pars->GC_11, amp[864]); 
  FFV1_0(w[45], w[117], w[4], pars->GC_11, amp[865]); 
  VVV1_0(w[96], w[6], w[144], pars->GC_10, amp[866]); 
  FFV1_0(w[3], w[117], w[96], pars->GC_11, amp[867]); 
  VVV1_0(w[4], w[46], w[144], pars->GC_10, amp[868]); 
  FFV1_0(w[3], w[130], w[46], pars->GC_11, amp[869]); 
  FFV1_0(w[3], w[116], w[163], pars->GC_11, amp[870]); 
  FFV1_0(w[3], w[116], w[164], pars->GC_11, amp[871]); 
  FFV1_0(w[3], w[116], w[165], pars->GC_11, amp[872]); 
  FFV1_0(w[57], w[127], w[6], pars->GC_11, amp[873]); 
  FFV1_0(w[12], w[127], w[4], pars->GC_11, amp[874]); 
  VVV1_0(w[96], w[6], w[148], pars->GC_10, amp[875]); 
  FFV1_0(w[12], w[2], w[96], pars->GC_11, amp[876]); 
  VVV1_0(w[4], w[46], w[148], pars->GC_10, amp[877]); 
  FFV1_0(w[57], w[2], w[46], pars->GC_11, amp[878]); 
  FFV1_0(w[10], w[2], w[163], pars->GC_11, amp[879]); 
  FFV1_0(w[10], w[2], w[164], pars->GC_11, amp[880]); 
  FFV1_0(w[10], w[2], w[165], pars->GC_11, amp[881]); 
  FFV1_0(w[45], w[131], w[6], pars->GC_11, amp[882]); 
  FFV1_0(w[45], w[121], w[4], pars->GC_11, amp[883]); 
  VVV1_0(w[96], w[6], w[149], pars->GC_10, amp[884]); 
  FFV1_0(w[3], w[121], w[96], pars->GC_11, amp[885]); 
  VVV1_0(w[4], w[46], w[149], pars->GC_10, amp[886]); 
  FFV1_0(w[3], w[131], w[46], pars->GC_11, amp[887]); 
  FFV1_0(w[3], w[120], w[163], pars->GC_11, amp[888]); 
  FFV1_0(w[3], w[120], w[164], pars->GC_11, amp[889]); 
  FFV1_0(w[3], w[120], w[165], pars->GC_11, amp[890]); 
  FFV1_0(w[59], w[127], w[6], pars->GC_11, amp[891]); 
  FFV1_0(w[19], w[127], w[4], pars->GC_11, amp[892]); 
  VVV1_0(w[96], w[6], w[150], pars->GC_10, amp[893]); 
  FFV1_0(w[19], w[2], w[96], pars->GC_11, amp[894]); 
  VVV1_0(w[4], w[46], w[150], pars->GC_10, amp[895]); 
  FFV1_0(w[59], w[2], w[46], pars->GC_11, amp[896]); 
  FFV1_0(w[18], w[2], w[163], pars->GC_11, amp[897]); 
  FFV1_0(w[18], w[2], w[164], pars->GC_11, amp[898]); 
  FFV1_0(w[18], w[2], w[165], pars->GC_11, amp[899]); 
  FFV1_0(w[48], w[130], w[5], pars->GC_11, amp[900]); 
  FFV1_0(w[48], w[119], w[4], pars->GC_11, amp[901]); 
  VVV1_0(w[77], w[5], w[144], pars->GC_10, amp[902]); 
  FFV1_0(w[3], w[119], w[77], pars->GC_11, amp[903]); 
  VVV1_0(w[4], w[49], w[144], pars->GC_10, amp[904]); 
  FFV1_0(w[3], w[130], w[49], pars->GC_11, amp[905]); 
  FFV1_0(w[3], w[116], w[166], pars->GC_11, amp[906]); 
  FFV1_0(w[3], w[116], w[167], pars->GC_11, amp[907]); 
  FFV1_0(w[3], w[116], w[168], pars->GC_11, amp[908]); 
  FFV1_0(w[57], w[129], w[5], pars->GC_11, amp[909]); 
  FFV1_0(w[15], w[129], w[4], pars->GC_11, amp[910]); 
  VVV1_0(w[77], w[5], w[148], pars->GC_10, amp[911]); 
  FFV1_0(w[15], w[2], w[77], pars->GC_11, amp[912]); 
  VVV1_0(w[4], w[49], w[148], pars->GC_10, amp[913]); 
  FFV1_0(w[57], w[2], w[49], pars->GC_11, amp[914]); 
  FFV1_0(w[10], w[2], w[166], pars->GC_11, amp[915]); 
  FFV1_0(w[10], w[2], w[167], pars->GC_11, amp[916]); 
  FFV1_0(w[10], w[2], w[168], pars->GC_11, amp[917]); 
  FFV1_0(w[48], w[131], w[5], pars->GC_11, amp[918]); 
  FFV1_0(w[48], w[123], w[4], pars->GC_11, amp[919]); 
  VVV1_0(w[77], w[5], w[149], pars->GC_10, amp[920]); 
  FFV1_0(w[3], w[123], w[77], pars->GC_11, amp[921]); 
  VVV1_0(w[4], w[49], w[149], pars->GC_10, amp[922]); 
  FFV1_0(w[3], w[131], w[49], pars->GC_11, amp[923]); 
  FFV1_0(w[3], w[120], w[166], pars->GC_11, amp[924]); 
  FFV1_0(w[3], w[120], w[167], pars->GC_11, amp[925]); 
  FFV1_0(w[3], w[120], w[168], pars->GC_11, amp[926]); 
  FFV1_0(w[59], w[129], w[5], pars->GC_11, amp[927]); 
  FFV1_0(w[21], w[129], w[4], pars->GC_11, amp[928]); 
  VVV1_0(w[77], w[5], w[150], pars->GC_10, amp[929]); 
  FFV1_0(w[21], w[2], w[77], pars->GC_11, amp[930]); 
  VVV1_0(w[4], w[49], w[150], pars->GC_10, amp[931]); 
  FFV1_0(w[59], w[2], w[49], pars->GC_11, amp[932]); 
  FFV1_0(w[18], w[2], w[166], pars->GC_11, amp[933]); 
  FFV1_0(w[18], w[2], w[167], pars->GC_11, amp[934]); 
  FFV1_0(w[18], w[2], w[168], pars->GC_11, amp[935]); 
  FFV1_0(w[169], w[116], w[7], pars->GC_11, amp[936]); 
  FFV1_0(w[170], w[116], w[7], pars->GC_11, amp[937]); 
  FFV1_0(w[171], w[116], w[7], pars->GC_11, amp[938]); 
  FFV1_0(w[3], w[116], w[172], pars->GC_11, amp[939]); 
  FFV1_0(w[3], w[116], w[173], pars->GC_11, amp[940]); 
  FFV1_0(w[3], w[116], w[174], pars->GC_11, amp[941]); 
  FFV1_0(w[10], w[175], w[7], pars->GC_11, amp[942]); 
  FFV1_0(w[10], w[176], w[7], pars->GC_11, amp[943]); 
  FFV1_0(w[10], w[177], w[7], pars->GC_11, amp[944]); 
  FFV1_0(w[10], w[2], w[172], pars->GC_11, amp[945]); 
  FFV1_0(w[10], w[2], w[173], pars->GC_11, amp[946]); 
  FFV1_0(w[10], w[2], w[174], pars->GC_11, amp[947]); 
  FFV1_0(w[169], w[120], w[7], pars->GC_11, amp[948]); 
  FFV1_0(w[170], w[120], w[7], pars->GC_11, amp[949]); 
  FFV1_0(w[171], w[120], w[7], pars->GC_11, amp[950]); 
  FFV1_0(w[3], w[120], w[172], pars->GC_11, amp[951]); 
  FFV1_0(w[3], w[120], w[173], pars->GC_11, amp[952]); 
  FFV1_0(w[3], w[120], w[174], pars->GC_11, amp[953]); 
  FFV1_0(w[18], w[175], w[7], pars->GC_11, amp[954]); 
  FFV1_0(w[18], w[176], w[7], pars->GC_11, amp[955]); 
  FFV1_0(w[18], w[177], w[7], pars->GC_11, amp[956]); 
  FFV1_0(w[18], w[2], w[172], pars->GC_11, amp[957]); 
  FFV1_0(w[18], w[2], w[173], pars->GC_11, amp[958]); 
  FFV1_0(w[18], w[2], w[174], pars->GC_11, amp[959]); 
  FFV1_0(w[178], w[116], w[6], pars->GC_11, amp[960]); 
  FFV1_0(w[179], w[116], w[6], pars->GC_11, amp[961]); 
  FFV1_0(w[180], w[116], w[6], pars->GC_11, amp[962]); 
  FFV1_0(w[3], w[116], w[181], pars->GC_11, amp[963]); 
  FFV1_0(w[3], w[116], w[182], pars->GC_11, amp[964]); 
  FFV1_0(w[3], w[116], w[183], pars->GC_11, amp[965]); 
  FFV1_0(w[10], w[184], w[6], pars->GC_11, amp[966]); 
  FFV1_0(w[10], w[185], w[6], pars->GC_11, amp[967]); 
  FFV1_0(w[10], w[186], w[6], pars->GC_11, amp[968]); 
  FFV1_0(w[10], w[2], w[181], pars->GC_11, amp[969]); 
  FFV1_0(w[10], w[2], w[182], pars->GC_11, amp[970]); 
  FFV1_0(w[10], w[2], w[183], pars->GC_11, amp[971]); 
  FFV1_0(w[178], w[120], w[6], pars->GC_11, amp[972]); 
  FFV1_0(w[179], w[120], w[6], pars->GC_11, amp[973]); 
  FFV1_0(w[180], w[120], w[6], pars->GC_11, amp[974]); 
  FFV1_0(w[3], w[120], w[181], pars->GC_11, amp[975]); 
  FFV1_0(w[3], w[120], w[182], pars->GC_11, amp[976]); 
  FFV1_0(w[3], w[120], w[183], pars->GC_11, amp[977]); 
  FFV1_0(w[18], w[184], w[6], pars->GC_11, amp[978]); 
  FFV1_0(w[18], w[185], w[6], pars->GC_11, amp[979]); 
  FFV1_0(w[18], w[186], w[6], pars->GC_11, amp[980]); 
  FFV1_0(w[18], w[2], w[181], pars->GC_11, amp[981]); 
  FFV1_0(w[18], w[2], w[182], pars->GC_11, amp[982]); 
  FFV1_0(w[18], w[2], w[183], pars->GC_11, amp[983]); 
  FFV1_0(w[187], w[116], w[5], pars->GC_11, amp[984]); 
  FFV1_0(w[188], w[116], w[5], pars->GC_11, amp[985]); 
  FFV1_0(w[189], w[116], w[5], pars->GC_11, amp[986]); 
  FFV1_0(w[3], w[116], w[190], pars->GC_11, amp[987]); 
  FFV1_0(w[3], w[116], w[191], pars->GC_11, amp[988]); 
  FFV1_0(w[3], w[116], w[192], pars->GC_11, amp[989]); 
  FFV1_0(w[10], w[193], w[5], pars->GC_11, amp[990]); 
  FFV1_0(w[10], w[194], w[5], pars->GC_11, amp[991]); 
  FFV1_0(w[10], w[195], w[5], pars->GC_11, amp[992]); 
  FFV1_0(w[10], w[2], w[190], pars->GC_11, amp[993]); 
  FFV1_0(w[10], w[2], w[191], pars->GC_11, amp[994]); 
  FFV1_0(w[10], w[2], w[192], pars->GC_11, amp[995]); 
  FFV1_0(w[187], w[120], w[5], pars->GC_11, amp[996]); 
  FFV1_0(w[188], w[120], w[5], pars->GC_11, amp[997]); 
  FFV1_0(w[189], w[120], w[5], pars->GC_11, amp[998]); 
  FFV1_0(w[3], w[120], w[190], pars->GC_11, amp[999]); 
  FFV1_0(w[3], w[120], w[191], pars->GC_11, amp[1000]); 
  FFV1_0(w[3], w[120], w[192], pars->GC_11, amp[1001]); 
  FFV1_0(w[18], w[193], w[5], pars->GC_11, amp[1002]); 
  FFV1_0(w[18], w[194], w[5], pars->GC_11, amp[1003]); 
  FFV1_0(w[18], w[195], w[5], pars->GC_11, amp[1004]); 
  FFV1_0(w[18], w[2], w[190], pars->GC_11, amp[1005]); 
  FFV1_0(w[18], w[2], w[191], pars->GC_11, amp[1006]); 
  FFV1_0(w[18], w[2], w[192], pars->GC_11, amp[1007]); 
  FFV1_0(w[196], w[116], w[4], pars->GC_11, amp[1008]); 
  FFV1_0(w[197], w[116], w[4], pars->GC_11, amp[1009]); 
  FFV1_0(w[198], w[116], w[4], pars->GC_11, amp[1010]); 
  FFV1_0(w[3], w[116], w[199], pars->GC_11, amp[1011]); 
  FFV1_0(w[3], w[116], w[200], pars->GC_11, amp[1012]); 
  FFV1_0(w[3], w[116], w[201], pars->GC_11, amp[1013]); 
  FFV1_0(w[10], w[202], w[4], pars->GC_11, amp[1014]); 
  FFV1_0(w[10], w[203], w[4], pars->GC_11, amp[1015]); 
  FFV1_0(w[10], w[204], w[4], pars->GC_11, amp[1016]); 
  FFV1_0(w[10], w[2], w[199], pars->GC_11, amp[1017]); 
  FFV1_0(w[10], w[2], w[200], pars->GC_11, amp[1018]); 
  FFV1_0(w[10], w[2], w[201], pars->GC_11, amp[1019]); 
  FFV1_0(w[196], w[120], w[4], pars->GC_11, amp[1020]); 
  FFV1_0(w[197], w[120], w[4], pars->GC_11, amp[1021]); 
  FFV1_0(w[198], w[120], w[4], pars->GC_11, amp[1022]); 
  FFV1_0(w[3], w[120], w[199], pars->GC_11, amp[1023]); 
  FFV1_0(w[3], w[120], w[200], pars->GC_11, amp[1024]); 
  FFV1_0(w[3], w[120], w[201], pars->GC_11, amp[1025]); 
  FFV1_0(w[18], w[202], w[4], pars->GC_11, amp[1026]); 
  FFV1_0(w[18], w[203], w[4], pars->GC_11, amp[1027]); 
  FFV1_0(w[18], w[204], w[4], pars->GC_11, amp[1028]); 
  FFV1_0(w[18], w[2], w[199], pars->GC_11, amp[1029]); 
  FFV1_0(w[18], w[2], w[200], pars->GC_11, amp[1030]); 
  FFV1_0(w[18], w[2], w[201], pars->GC_11, amp[1031]); 

}
double eeuugggg::matrix_1_epem_uuxgggg()
{
//  int i, j;
  // Local variables
//  const int ngraphs = 1032; 
  const int ncolor = 24;
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  //static const double denom[ncolor] = {54, 54, 54, 54, 54, 54, 54, 54, 54, 54,
  //    54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54};
  /*static const double cf[ncolor][ncolor] = {{512, -64, -64, 8, 8, 80, -64, 8,
      8, -1, -1, -10, 8, -1, 80, -10, 71, 62, -1, -10, -10, 62, 62, -28}, {-64,
      512, 8, 80, -64, 8, 8, -64, -1, -10, 8, -1, -1, -10, -10, 62, 62, -28, 8,
      -1, 80, -10, 71, 62}, {-64, 8, 512, -64, 80, 8, 8, -1, 80, -10, 71, 62,
      -64, 8, 8, -1, -1, -10, -10, -1, 62, -28, -10, 62}, {8, 80, -64, 512, 8,
      -64, -1, -10, -10, 62, 62, -28, 8, -64, -1, -10, 8, -1, -1, 8, 71, 62,
      80, -10}, {8, -64, 80, 8, 512, -64, -1, 8, 71, 62, 80, -10, -10, -1, 62,
      -28, -10, 62, -64, 8, 8, -1, -1, -10}, {80, 8, 8, -64, -64, 512, -10, -1,
      62, -28, -10, 62, -1, 8, 71, 62, 80, -10, 8, -64, -1, -10, 8, -1}, {-64,
      8, 8, -1, -1, -10, 512, -64, -64, 8, 8, 80, 80, -10, 8, -1, 62, 71, -10,
      62, -1, -10, -28, 62}, {8, -64, -1, -10, 8, -1, -64, 512, 8, 80, -64, 8,
      -10, 62, -1, -10, -28, 62, 80, -10, 8, -1, 62, 71}, {8, -1, 80, -10, 71,
      62, -64, 8, 512, -64, 80, 8, 8, -1, -64, 8, -10, -1, 62, -28, -10, -1,
      62, -10}, {-1, -10, -10, 62, 62, -28, 8, 80, -64, 512, 8, -64, -1, -10,
      8, -64, -1, 8, 71, 62, -1, 8, -10, 80}, {-1, 8, 71, 62, 80, -10, 8, -64,
      80, 8, 512, -64, 62, -28, -10, -1, 62, -10, 8, -1, -64, 8, -10, -1},
      {-10, -1, 62, -28, -10, 62, 80, 8, 8, -64, -64, 512, 71, 62, -1, 8, -10,
      80, -1, -10, 8, -64, -1, 8}, {8, -1, -64, 8, -10, -1, 80, -10, 8, -1, 62,
      71, 512, -64, -64, 8, 8, 80, 62, -10, -28, 62, -1, -10}, {-1, -10, 8,
      -64, -1, 8, -10, 62, -1, -10, -28, 62, -64, 512, 8, 80, -64, 8, -10, 80,
      62, 71, 8, -1}, {80, -10, 8, -1, 62, 71, 8, -1, -64, 8, -10, -1, -64, 8,
      512, -64, 80, 8, -28, 62, 62, -10, -10, -1}, {-10, 62, -1, -10, -28, 62,
      -1, -10, 8, -64, -1, 8, 8, 80, -64, 512, 8, -64, 62, 71, -10, 80, -1, 8},
      {71, 62, -1, 8, -10, 80, 62, -28, -10, -1, 62, -10, 8, -64, 80, 8, 512,
      -64, -1, 8, -10, -1, -64, 8}, {62, -28, -10, -1, 62, -10, 71, 62, -1, 8,
      -10, 80, 80, 8, 8, -64, -64, 512, -10, -1, -1, 8, 8, -64}, {-1, 8, -10,
      -1, -64, 8, -10, 80, 62, 71, 8, -1, 62, -10, -28, 62, -1, -10, 512, -64,
      -64, 8, 8, 80}, {-10, -1, -1, 8, 8, -64, 62, -10, -28, 62, -1, -10, -10,
      80, 62, 71, 8, -1, -64, 512, 8, 80, -64, 8}, {-10, 80, 62, 71, 8, -1, -1,
      8, -10, -1, -64, 8, -28, 62, 62, -10, -10, -1, -64, 8, 512, -64, 80, 8},
      {62, -10, -28, 62, -1, -10, -10, -1, -1, 8, 8, -64, 62, 71, -10, 80, -1,
      8, 8, 80, -64, 512, 8, -64}, {62, 71, -10, 80, -1, 8, -28, 62, 62, -10,
      -10, -1, -1, 8, -10, -1, -64, 8, 8, -64, 80, 8, 512, -64}, {-28, 62, 62,
      -10, -10, -1, 62, 71, -10, 80, -1, 8, -10, -1, -1, 8, 8, -64, 80, 8, 8,
      -64, -64, 512}};
   */

  // Calculate color flows
  jamp[0] = +amp[1] + amp[7] + amp[45] + amp[46] + amp[48] + amp[51] + amp[52]
      + amp[54] - std::complex<double> (0, 1) * amp[56] - std::complex<double>
      (0, 1) * amp[57] - std::complex<double> (0, 1) * amp[58] -
      std::complex<double> (0, 1) * amp[59] - amp[61] - std::complex<double>
      (0, 1) * amp[62] - amp[64] - amp[67] - std::complex<double> (0, 1) *
      amp[68] - amp[70] - std::complex<double> (0, 1) * amp[84] - amp[85] -
      std::complex<double> (0, 1) * amp[86] - amp[88] - std::complex<double>
      (0, 1) * amp[89] - std::complex<double> (0, 1) * amp[90] - amp[91] -
      std::complex<double> (0, 1) * amp[92] - amp[94] - std::complex<double>
      (0, 1) * amp[95] + amp[98] - amp[96] + amp[101] - amp[99] + amp[104] -
      amp[102] + amp[107] - amp[105] + amp[616] + amp[622] -
      std::complex<double> (0, 1) * amp[625] - amp[626] - std::complex<double>
      (0, 1) * amp[627] - amp[628] - std::complex<double> (0, 1) * amp[629] -
      std::complex<double> (0, 1) * amp[631] - amp[632] - std::complex<double>
      (0, 1) * amp[633] - amp[634] - std::complex<double> (0, 1) * amp[635] -
      std::complex<double> (0, 1) * amp[649] - amp[650] - amp[652] -
      std::complex<double> (0, 1) * amp[655] - amp[656] - amp[658] + amp[662] -
      amp[660] + amp[665] - amp[663] + amp[668] - amp[666] + amp[671] -
      amp[669] + std::complex<double> (0, 1) * amp[674] - std::complex<double>
      (0, 1) * amp[680] + std::complex<double> (0, 1) * amp[678] -
      std::complex<double> (0, 1) * amp[682] + std::complex<double> (0, 1) *
      amp[683] - amp[684] - std::complex<double> (0, 1) * amp[689] +
      std::complex<double> (0, 1) * amp[687] + std::complex<double> (0, 1) *
      amp[692] - std::complex<double> (0, 1) * amp[698] + std::complex<double>
      (0, 1) * amp[696] - std::complex<double> (0, 1) * amp[700] +
      std::complex<double> (0, 1) * amp[701] - amp[702] - std::complex<double>
      (0, 1) * amp[707] + std::complex<double> (0, 1) * amp[705] - amp[709] +
      std::complex<double> (0, 1) * amp[710] - amp[711] + std::complex<double>
      (0, 1) * amp[713] - amp[714] - amp[717] + std::complex<double> (0, 1) *
      amp[718] - amp[719] + std::complex<double> (0, 1) * amp[721] - amp[722] +
      std::complex<double> (0, 1) * amp[830] + std::complex<double> (0, 1) *
      amp[832] - amp[833] - std::complex<double> (0, 1) * amp[836] +
      std::complex<double> (0, 1) * amp[834] + std::complex<double> (0, 1) *
      amp[839] - amp[840] + std::complex<double> (0, 1) * amp[841] -
      std::complex<double> (0, 1) * amp[845] + std::complex<double> (0, 1) *
      amp[843] + std::complex<double> (0, 1) * amp[848] + std::complex<double>
      (0, 1) * amp[850] - amp[851] - std::complex<double> (0, 1) * amp[854] +
      std::complex<double> (0, 1) * amp[852] + std::complex<double> (0, 1) *
      amp[857] - amp[858] + std::complex<double> (0, 1) * amp[859] -
      std::complex<double> (0, 1) * amp[863] + std::complex<double> (0, 1) *
      amp[861] - std::complex<double> (0, 1) * amp[900] + std::complex<double>
      (0, 1) * amp[904] - amp[905] - std::complex<double> (0, 1) * amp[908] +
      std::complex<double> (0, 1) * amp[906] + std::complex<double> (0, 1) *
      amp[913] - std::complex<double> (0, 1) * amp[917] + std::complex<double>
      (0, 1) * amp[915] - std::complex<double> (0, 1) * amp[918] +
      std::complex<double> (0, 1) * amp[922] - amp[923] - std::complex<double>
      (0, 1) * amp[926] + std::complex<double> (0, 1) * amp[924] +
      std::complex<double> (0, 1) * amp[931] - std::complex<double> (0, 1) *
      amp[935] + std::complex<double> (0, 1) * amp[933] - std::complex<double>
      (0, 1) * amp[941] + std::complex<double> (0, 1) * amp[939] + amp[944] -
      amp[942] - std::complex<double> (0, 1) * amp[947] + std::complex<double>
      (0, 1) * amp[945] - std::complex<double> (0, 1) * amp[953] +
      std::complex<double> (0, 1) * amp[951] + amp[956] - amp[954] -
      std::complex<double> (0, 1) * amp[959] + std::complex<double> (0, 1) *
      amp[957] + amp[1010] - amp[1008] - std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1011] -
      std::complex<double> (0, 1) * amp[1019] + std::complex<double> (0, 1) *
      amp[1017] + amp[1022] - amp[1020] - std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1023] -
      std::complex<double> (0, 1) * amp[1031] + std::complex<double> (0, 1) *
      amp[1029];
  jamp[1] = +amp[0] + amp[6] + amp[29] + amp[30] + amp[32] + amp[35] + amp[36]
      + amp[38] - std::complex<double> (0, 1) * amp[40] - std::complex<double>
      (0, 1) * amp[41] - std::complex<double> (0, 1) * amp[42] -
      std::complex<double> (0, 1) * amp[43] - amp[73] - std::complex<double>
      (0, 1) * amp[74] - amp[76] - amp[79] - std::complex<double> (0, 1) *
      amp[80] - amp[82] + std::complex<double> (0, 1) * amp[84] + amp[85] +
      std::complex<double> (0, 1) * amp[86] + amp[88] + std::complex<double>
      (0, 1) * amp[89] + std::complex<double> (0, 1) * amp[90] + amp[91] +
      std::complex<double> (0, 1) * amp[92] + amp[94] + std::complex<double>
      (0, 1) * amp[95] + amp[96] + amp[97] + amp[99] + amp[100] + amp[102] +
      amp[103] + amp[105] + amp[106] + amp[556] + amp[562] -
      std::complex<double> (0, 1) * amp[565] - amp[566] - std::complex<double>
      (0, 1) * amp[567] - amp[568] - std::complex<double> (0, 1) * amp[569] -
      std::complex<double> (0, 1) * amp[571] - amp[572] - std::complex<double>
      (0, 1) * amp[573] - amp[574] - std::complex<double> (0, 1) * amp[575] -
      std::complex<double> (0, 1) * amp[589] - amp[590] - amp[592] -
      std::complex<double> (0, 1) * amp[595] - amp[596] - amp[598] + amp[602] -
      amp[600] + amp[605] - amp[603] + amp[608] - amp[606] + amp[611] -
      amp[609] + std::complex<double> (0, 1) * amp[676] - std::complex<double>
      (0, 1) * amp[678] - std::complex<double> (0, 1) * amp[679] -
      std::complex<double> (0, 1) * amp[681] + std::complex<double> (0, 1) *
      amp[685] - amp[686] - std::complex<double> (0, 1) * amp[687] -
      std::complex<double> (0, 1) * amp[688] + std::complex<double> (0, 1) *
      amp[694] - std::complex<double> (0, 1) * amp[696] - std::complex<double>
      (0, 1) * amp[697] - std::complex<double> (0, 1) * amp[699] +
      std::complex<double> (0, 1) * amp[703] - amp[704] - std::complex<double>
      (0, 1) * amp[705] - std::complex<double> (0, 1) * amp[706] + amp[709] -
      std::complex<double> (0, 1) * amp[710] + amp[711] - std::complex<double>
      (0, 1) * amp[713] + amp[714] + amp[717] - std::complex<double> (0, 1) *
      amp[718] + amp[719] - std::complex<double> (0, 1) * amp[721] + amp[722] +
      std::complex<double> (0, 1) * amp[866] + std::complex<double> (0, 1) *
      amp[868] - amp[869] - std::complex<double> (0, 1) * amp[872] +
      std::complex<double> (0, 1) * amp[870] + std::complex<double> (0, 1) *
      amp[875] - amp[876] + std::complex<double> (0, 1) * amp[877] -
      std::complex<double> (0, 1) * amp[881] + std::complex<double> (0, 1) *
      amp[879] + std::complex<double> (0, 1) * amp[884] + std::complex<double>
      (0, 1) * amp[886] - amp[887] - std::complex<double> (0, 1) * amp[890] +
      std::complex<double> (0, 1) * amp[888] + std::complex<double> (0, 1) *
      amp[893] - amp[894] + std::complex<double> (0, 1) * amp[895] -
      std::complex<double> (0, 1) * amp[899] + std::complex<double> (0, 1) *
      amp[897] + std::complex<double> (0, 1) * amp[900] - std::complex<double>
      (0, 1) * amp[904] + amp[905] + std::complex<double> (0, 1) * amp[908] -
      std::complex<double> (0, 1) * amp[906] - std::complex<double> (0, 1) *
      amp[913] + std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[915] + std::complex<double> (0, 1) * amp[918] -
      std::complex<double> (0, 1) * amp[922] + amp[923] + std::complex<double>
      (0, 1) * amp[926] - std::complex<double> (0, 1) * amp[924] -
      std::complex<double> (0, 1) * amp[931] + std::complex<double> (0, 1) *
      amp[935] - std::complex<double> (0, 1) * amp[933] - std::complex<double>
      (0, 1) * amp[965] + std::complex<double> (0, 1) * amp[963] + amp[968] -
      amp[966] - std::complex<double> (0, 1) * amp[971] + std::complex<double>
      (0, 1) * amp[969] - std::complex<double> (0, 1) * amp[977] +
      std::complex<double> (0, 1) * amp[975] + amp[980] - amp[978] -
      std::complex<double> (0, 1) * amp[983] + std::complex<double> (0, 1) *
      amp[981] + amp[1008] + amp[1009] - std::complex<double> (0, 1) *
      amp[1011] - std::complex<double> (0, 1) * amp[1012] -
      std::complex<double> (0, 1) * amp[1017] - std::complex<double> (0, 1) *
      amp[1018] + amp[1020] + amp[1021] - std::complex<double> (0, 1) *
      amp[1023] - std::complex<double> (0, 1) * amp[1024] -
      std::complex<double> (0, 1) * amp[1029] - std::complex<double> (0, 1) *
      amp[1030];
  jamp[2] = +amp[3] + amp[9] + amp[44] + amp[47] + amp[49] + amp[50] + amp[53]
      + amp[55] + std::complex<double> (0, 1) * amp[56] + std::complex<double>
      (0, 1) * amp[57] + std::complex<double> (0, 1) * amp[58] +
      std::complex<double> (0, 1) * amp[59] + amp[61] + std::complex<double>
      (0, 1) * amp[62] + amp[64] + amp[67] + std::complex<double> (0, 1) *
      amp[68] + amp[70] - std::complex<double> (0, 1) * amp[72] + amp[73] -
      std::complex<double> (0, 1) * amp[75] + amp[76] - std::complex<double>
      (0, 1) * amp[77] - std::complex<double> (0, 1) * amp[78] + amp[79] -
      std::complex<double> (0, 1) * amp[81] + amp[82] - std::complex<double>
      (0, 1) * amp[83] - amp[98] - amp[97] - amp[101] - amp[100] - amp[104] -
      amp[103] - amp[107] - amp[106] + amp[614] + amp[620] -
      std::complex<double> (0, 1) * amp[637] - amp[638] - std::complex<double>
      (0, 1) * amp[639] - amp[640] - std::complex<double> (0, 1) * amp[641] -
      std::complex<double> (0, 1) * amp[643] - amp[644] - std::complex<double>
      (0, 1) * amp[645] - amp[646] - std::complex<double> (0, 1) * amp[647] +
      std::complex<double> (0, 1) * amp[649] + amp[650] + amp[652] +
      std::complex<double> (0, 1) * amp[655] + amp[656] + amp[658] + amp[660] +
      amp[661] + amp[663] + amp[664] + amp[666] + amp[667] + amp[669] +
      amp[670] + std::complex<double> (0, 1) * amp[726] - std::complex<double>
      (0, 1) * amp[732] + std::complex<double> (0, 1) * amp[730] -
      std::complex<double> (0, 1) * amp[734] + std::complex<double> (0, 1) *
      amp[735] - amp[736] - std::complex<double> (0, 1) * amp[741] +
      std::complex<double> (0, 1) * amp[739] + std::complex<double> (0, 1) *
      amp[744] - std::complex<double> (0, 1) * amp[750] + std::complex<double>
      (0, 1) * amp[748] - std::complex<double> (0, 1) * amp[752] +
      std::complex<double> (0, 1) * amp[753] - amp[754] - std::complex<double>
      (0, 1) * amp[759] + std::complex<double> (0, 1) * amp[757] - amp[761] +
      std::complex<double> (0, 1) * amp[762] - amp[763] + std::complex<double>
      (0, 1) * amp[765] - amp[766] - amp[769] + std::complex<double> (0, 1) *
      amp[770] - amp[771] + std::complex<double> (0, 1) * amp[773] - amp[774] -
      std::complex<double> (0, 1) * amp[830] - std::complex<double> (0, 1) *
      amp[832] + amp[833] + std::complex<double> (0, 1) * amp[836] -
      std::complex<double> (0, 1) * amp[834] - std::complex<double> (0, 1) *
      amp[839] + amp[840] - std::complex<double> (0, 1) * amp[841] +
      std::complex<double> (0, 1) * amp[845] - std::complex<double> (0, 1) *
      amp[843] - std::complex<double> (0, 1) * amp[848] - std::complex<double>
      (0, 1) * amp[850] + amp[851] + std::complex<double> (0, 1) * amp[854] -
      std::complex<double> (0, 1) * amp[852] - std::complex<double> (0, 1) *
      amp[857] + amp[858] - std::complex<double> (0, 1) * amp[859] +
      std::complex<double> (0, 1) * amp[863] - std::complex<double> (0, 1) *
      amp[861] - std::complex<double> (0, 1) * amp[864] - std::complex<double>
      (0, 1) * amp[868] + amp[869] - std::complex<double> (0, 1) * amp[870] -
      std::complex<double> (0, 1) * amp[871] - std::complex<double> (0, 1) *
      amp[877] - std::complex<double> (0, 1) * amp[879] - std::complex<double>
      (0, 1) * amp[880] - std::complex<double> (0, 1) * amp[882] -
      std::complex<double> (0, 1) * amp[886] + amp[887] - std::complex<double>
      (0, 1) * amp[888] - std::complex<double> (0, 1) * amp[889] -
      std::complex<double> (0, 1) * amp[895] - std::complex<double> (0, 1) *
      amp[897] - std::complex<double> (0, 1) * amp[898] - std::complex<double>
      (0, 1) * amp[939] - std::complex<double> (0, 1) * amp[940] + amp[942] +
      amp[943] - std::complex<double> (0, 1) * amp[945] - std::complex<double>
      (0, 1) * amp[946] - std::complex<double> (0, 1) * amp[951] -
      std::complex<double> (0, 1) * amp[952] + amp[954] + amp[955] -
      std::complex<double> (0, 1) * amp[957] - std::complex<double> (0, 1) *
      amp[958] - amp[1010] - amp[1009] + std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1012] +
      std::complex<double> (0, 1) * amp[1019] + std::complex<double> (0, 1) *
      amp[1018] - amp[1022] - amp[1021] + std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1024] +
      std::complex<double> (0, 1) * amp[1031] + std::complex<double> (0, 1) *
      amp[1030];
  jamp[3] = +amp[2] + amp[8] + amp[13] + amp[14] + amp[16] + amp[19] + amp[20]
      + amp[22] - std::complex<double> (0, 1) * amp[24] - std::complex<double>
      (0, 1) * amp[25] - std::complex<double> (0, 1) * amp[26] -
      std::complex<double> (0, 1) * amp[27] + std::complex<double> (0, 1) *
      amp[72] - amp[73] + std::complex<double> (0, 1) * amp[75] - amp[76] +
      std::complex<double> (0, 1) * amp[77] + std::complex<double> (0, 1) *
      amp[78] - amp[79] + std::complex<double> (0, 1) * amp[81] - amp[82] +
      std::complex<double> (0, 1) * amp[83] + amp[85] - std::complex<double>
      (0, 1) * amp[87] + amp[88] + amp[91] - std::complex<double> (0, 1) *
      amp[93] + amp[94] + amp[96] + amp[97] + amp[99] + amp[100] + amp[102] +
      amp[103] + amp[105] + amp[106] + amp[496] + amp[502] -
      std::complex<double> (0, 1) * amp[505] - amp[506] - std::complex<double>
      (0, 1) * amp[507] - amp[508] - std::complex<double> (0, 1) * amp[509] -
      std::complex<double> (0, 1) * amp[511] - amp[512] - std::complex<double>
      (0, 1) * amp[513] - amp[514] - std::complex<double> (0, 1) * amp[515] -
      std::complex<double> (0, 1) * amp[529] - amp[530] - amp[532] -
      std::complex<double> (0, 1) * amp[535] - amp[536] - amp[538] + amp[542] -
      amp[540] + amp[545] - amp[543] + amp[548] - amp[546] + amp[551] -
      amp[549] + std::complex<double> (0, 1) * amp[728] - std::complex<double>
      (0, 1) * amp[730] - std::complex<double> (0, 1) * amp[731] -
      std::complex<double> (0, 1) * amp[733] + std::complex<double> (0, 1) *
      amp[737] - amp[738] - std::complex<double> (0, 1) * amp[739] -
      std::complex<double> (0, 1) * amp[740] + std::complex<double> (0, 1) *
      amp[746] - std::complex<double> (0, 1) * amp[748] - std::complex<double>
      (0, 1) * amp[749] - std::complex<double> (0, 1) * amp[751] +
      std::complex<double> (0, 1) * amp[755] - amp[756] - std::complex<double>
      (0, 1) * amp[757] - std::complex<double> (0, 1) * amp[758] + amp[761] -
      std::complex<double> (0, 1) * amp[762] + amp[763] - std::complex<double>
      (0, 1) * amp[765] + amp[766] + amp[769] - std::complex<double> (0, 1) *
      amp[770] + amp[771] - std::complex<double> (0, 1) * amp[773] + amp[774] +
      std::complex<double> (0, 1) * amp[864] + std::complex<double> (0, 1) *
      amp[868] - amp[869] + std::complex<double> (0, 1) * amp[870] +
      std::complex<double> (0, 1) * amp[871] + std::complex<double> (0, 1) *
      amp[877] + std::complex<double> (0, 1) * amp[879] + std::complex<double>
      (0, 1) * amp[880] + std::complex<double> (0, 1) * amp[882] +
      std::complex<double> (0, 1) * amp[886] - amp[887] + std::complex<double>
      (0, 1) * amp[888] + std::complex<double> (0, 1) * amp[889] +
      std::complex<double> (0, 1) * amp[895] + std::complex<double> (0, 1) *
      amp[897] + std::complex<double> (0, 1) * amp[898] + std::complex<double>
      (0, 1) * amp[902] - std::complex<double> (0, 1) * amp[904] + amp[905] -
      std::complex<double> (0, 1) * amp[906] - std::complex<double> (0, 1) *
      amp[907] + std::complex<double> (0, 1) * amp[911] - amp[912] -
      std::complex<double> (0, 1) * amp[913] - std::complex<double> (0, 1) *
      amp[915] - std::complex<double> (0, 1) * amp[916] + std::complex<double>
      (0, 1) * amp[920] - std::complex<double> (0, 1) * amp[922] + amp[923] -
      std::complex<double> (0, 1) * amp[924] - std::complex<double> (0, 1) *
      amp[925] + std::complex<double> (0, 1) * amp[929] - amp[930] -
      std::complex<double> (0, 1) * amp[931] - std::complex<double> (0, 1) *
      amp[933] - std::complex<double> (0, 1) * amp[934] - std::complex<double>
      (0, 1) * amp[989] + std::complex<double> (0, 1) * amp[987] + amp[992] -
      amp[990] - std::complex<double> (0, 1) * amp[995] + std::complex<double>
      (0, 1) * amp[993] - std::complex<double> (0, 1) * amp[1001] +
      std::complex<double> (0, 1) * amp[999] + amp[1004] - amp[1002] -
      std::complex<double> (0, 1) * amp[1007] + std::complex<double> (0, 1) *
      amp[1005] + amp[1008] + amp[1009] - std::complex<double> (0, 1) *
      amp[1011] - std::complex<double> (0, 1) * amp[1012] -
      std::complex<double> (0, 1) * amp[1017] - std::complex<double> (0, 1) *
      amp[1018] + amp[1020] + amp[1021] - std::complex<double> (0, 1) *
      amp[1023] - std::complex<double> (0, 1) * amp[1024] -
      std::complex<double> (0, 1) * amp[1029] - std::complex<double> (0, 1) *
      amp[1030];
  jamp[4] = +amp[5] + amp[11] + amp[28] + amp[31] + amp[33] + amp[34] + amp[37]
      + amp[39] + std::complex<double> (0, 1) * amp[40] + std::complex<double>
      (0, 1) * amp[41] + std::complex<double> (0, 1) * amp[42] +
      std::complex<double> (0, 1) * amp[43] - std::complex<double> (0, 1) *
      amp[60] + amp[61] - std::complex<double> (0, 1) * amp[63] + amp[64] -
      std::complex<double> (0, 1) * amp[65] - std::complex<double> (0, 1) *
      amp[66] + amp[67] - std::complex<double> (0, 1) * amp[69] + amp[70] -
      std::complex<double> (0, 1) * amp[71] + amp[73] + std::complex<double>
      (0, 1) * amp[74] + amp[76] + amp[79] + std::complex<double> (0, 1) *
      amp[80] + amp[82] - amp[98] - amp[97] - amp[101] - amp[100] - amp[104] -
      amp[103] - amp[107] - amp[106] + amp[554] + amp[560] -
      std::complex<double> (0, 1) * amp[577] - amp[578] - std::complex<double>
      (0, 1) * amp[579] - amp[580] - std::complex<double> (0, 1) * amp[581] -
      std::complex<double> (0, 1) * amp[583] - amp[584] - std::complex<double>
      (0, 1) * amp[585] - amp[586] - std::complex<double> (0, 1) * amp[587] +
      std::complex<double> (0, 1) * amp[589] + amp[590] + amp[592] +
      std::complex<double> (0, 1) * amp[595] + amp[596] + amp[598] + amp[600] +
      amp[601] + amp[603] + amp[604] + amp[606] + amp[607] + amp[609] +
      amp[610] + std::complex<double> (0, 1) * amp[778] - std::complex<double>
      (0, 1) * amp[784] + std::complex<double> (0, 1) * amp[782] -
      std::complex<double> (0, 1) * amp[786] + std::complex<double> (0, 1) *
      amp[787] - amp[788] - std::complex<double> (0, 1) * amp[793] +
      std::complex<double> (0, 1) * amp[791] + std::complex<double> (0, 1) *
      amp[796] - std::complex<double> (0, 1) * amp[802] + std::complex<double>
      (0, 1) * amp[800] - std::complex<double> (0, 1) * amp[804] +
      std::complex<double> (0, 1) * amp[805] - amp[806] - std::complex<double>
      (0, 1) * amp[811] + std::complex<double> (0, 1) * amp[809] - amp[813] +
      std::complex<double> (0, 1) * amp[814] - amp[815] + std::complex<double>
      (0, 1) * amp[817] - amp[818] - amp[821] + std::complex<double> (0, 1) *
      amp[822] - amp[823] + std::complex<double> (0, 1) * amp[825] - amp[826] -
      std::complex<double> (0, 1) * amp[828] - std::complex<double> (0, 1) *
      amp[832] + amp[833] - std::complex<double> (0, 1) * amp[834] -
      std::complex<double> (0, 1) * amp[835] - std::complex<double> (0, 1) *
      amp[841] - std::complex<double> (0, 1) * amp[843] - std::complex<double>
      (0, 1) * amp[844] - std::complex<double> (0, 1) * amp[846] -
      std::complex<double> (0, 1) * amp[850] + amp[851] - std::complex<double>
      (0, 1) * amp[852] - std::complex<double> (0, 1) * amp[853] -
      std::complex<double> (0, 1) * amp[859] - std::complex<double> (0, 1) *
      amp[861] - std::complex<double> (0, 1) * amp[862] - std::complex<double>
      (0, 1) * amp[866] - std::complex<double> (0, 1) * amp[868] + amp[869] +
      std::complex<double> (0, 1) * amp[872] - std::complex<double> (0, 1) *
      amp[870] - std::complex<double> (0, 1) * amp[875] + amp[876] -
      std::complex<double> (0, 1) * amp[877] + std::complex<double> (0, 1) *
      amp[881] - std::complex<double> (0, 1) * amp[879] - std::complex<double>
      (0, 1) * amp[884] - std::complex<double> (0, 1) * amp[886] + amp[887] +
      std::complex<double> (0, 1) * amp[890] - std::complex<double> (0, 1) *
      amp[888] - std::complex<double> (0, 1) * amp[893] + amp[894] -
      std::complex<double> (0, 1) * amp[895] + std::complex<double> (0, 1) *
      amp[899] - std::complex<double> (0, 1) * amp[897] - std::complex<double>
      (0, 1) * amp[963] - std::complex<double> (0, 1) * amp[964] + amp[966] +
      amp[967] - std::complex<double> (0, 1) * amp[969] - std::complex<double>
      (0, 1) * amp[970] - std::complex<double> (0, 1) * amp[975] -
      std::complex<double> (0, 1) * amp[976] + amp[978] + amp[979] -
      std::complex<double> (0, 1) * amp[981] - std::complex<double> (0, 1) *
      amp[982] - amp[1010] - amp[1009] + std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1012] +
      std::complex<double> (0, 1) * amp[1019] + std::complex<double> (0, 1) *
      amp[1018] - amp[1022] - amp[1021] + std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1024] +
      std::complex<double> (0, 1) * amp[1031] + std::complex<double> (0, 1) *
      amp[1030];
  jamp[5] = +amp[4] + amp[10] + amp[12] + amp[15] + amp[17] + amp[18] + amp[21]
      + amp[23] + std::complex<double> (0, 1) * amp[24] + std::complex<double>
      (0, 1) * amp[25] + std::complex<double> (0, 1) * amp[26] +
      std::complex<double> (0, 1) * amp[27] + std::complex<double> (0, 1) *
      amp[60] - amp[61] + std::complex<double> (0, 1) * amp[63] - amp[64] +
      std::complex<double> (0, 1) * amp[65] + std::complex<double> (0, 1) *
      amp[66] - amp[67] + std::complex<double> (0, 1) * amp[69] - amp[70] +
      std::complex<double> (0, 1) * amp[71] - amp[85] + std::complex<double>
      (0, 1) * amp[87] - amp[88] - amp[91] + std::complex<double> (0, 1) *
      amp[93] - amp[94] + amp[98] - amp[96] + amp[101] - amp[99] + amp[104] -
      amp[102] + amp[107] - amp[105] + amp[494] + amp[500] -
      std::complex<double> (0, 1) * amp[517] - amp[518] - std::complex<double>
      (0, 1) * amp[519] - amp[520] - std::complex<double> (0, 1) * amp[521] -
      std::complex<double> (0, 1) * amp[523] - amp[524] - std::complex<double>
      (0, 1) * amp[525] - amp[526] - std::complex<double> (0, 1) * amp[527] +
      std::complex<double> (0, 1) * amp[529] + amp[530] + amp[532] +
      std::complex<double> (0, 1) * amp[535] + amp[536] + amp[538] + amp[540] +
      amp[541] + amp[543] + amp[544] + amp[546] + amp[547] + amp[549] +
      amp[550] + std::complex<double> (0, 1) * amp[780] - std::complex<double>
      (0, 1) * amp[782] - std::complex<double> (0, 1) * amp[783] -
      std::complex<double> (0, 1) * amp[785] + std::complex<double> (0, 1) *
      amp[789] - amp[790] - std::complex<double> (0, 1) * amp[791] -
      std::complex<double> (0, 1) * amp[792] + std::complex<double> (0, 1) *
      amp[798] - std::complex<double> (0, 1) * amp[800] - std::complex<double>
      (0, 1) * amp[801] - std::complex<double> (0, 1) * amp[803] +
      std::complex<double> (0, 1) * amp[807] - amp[808] - std::complex<double>
      (0, 1) * amp[809] - std::complex<double> (0, 1) * amp[810] + amp[813] -
      std::complex<double> (0, 1) * amp[814] + amp[815] - std::complex<double>
      (0, 1) * amp[817] + amp[818] + amp[821] - std::complex<double> (0, 1) *
      amp[822] + amp[823] - std::complex<double> (0, 1) * amp[825] + amp[826] +
      std::complex<double> (0, 1) * amp[828] + std::complex<double> (0, 1) *
      amp[832] - amp[833] + std::complex<double> (0, 1) * amp[834] +
      std::complex<double> (0, 1) * amp[835] + std::complex<double> (0, 1) *
      amp[841] + std::complex<double> (0, 1) * amp[843] + std::complex<double>
      (0, 1) * amp[844] + std::complex<double> (0, 1) * amp[846] +
      std::complex<double> (0, 1) * amp[850] - amp[851] + std::complex<double>
      (0, 1) * amp[852] + std::complex<double> (0, 1) * amp[853] +
      std::complex<double> (0, 1) * amp[859] + std::complex<double> (0, 1) *
      amp[861] + std::complex<double> (0, 1) * amp[862] - std::complex<double>
      (0, 1) * amp[902] + std::complex<double> (0, 1) * amp[904] - amp[905] +
      std::complex<double> (0, 1) * amp[906] + std::complex<double> (0, 1) *
      amp[907] - std::complex<double> (0, 1) * amp[911] + amp[912] +
      std::complex<double> (0, 1) * amp[913] + std::complex<double> (0, 1) *
      amp[915] + std::complex<double> (0, 1) * amp[916] - std::complex<double>
      (0, 1) * amp[920] + std::complex<double> (0, 1) * amp[922] - amp[923] +
      std::complex<double> (0, 1) * amp[924] + std::complex<double> (0, 1) *
      amp[925] - std::complex<double> (0, 1) * amp[929] + amp[930] +
      std::complex<double> (0, 1) * amp[931] + std::complex<double> (0, 1) *
      amp[933] + std::complex<double> (0, 1) * amp[934] - std::complex<double>
      (0, 1) * amp[987] - std::complex<double> (0, 1) * amp[988] + amp[990] +
      amp[991] - std::complex<double> (0, 1) * amp[993] - std::complex<double>
      (0, 1) * amp[994] - std::complex<double> (0, 1) * amp[999] -
      std::complex<double> (0, 1) * amp[1000] + amp[1002] + amp[1003] -
      std::complex<double> (0, 1) * amp[1005] - std::complex<double> (0, 1) *
      amp[1006] + amp[1010] - amp[1008] - std::complex<double> (0, 1) *
      amp[1013] + std::complex<double> (0, 1) * amp[1011] -
      std::complex<double> (0, 1) * amp[1019] + std::complex<double> (0, 1) *
      amp[1017] + amp[1022] - amp[1020] - std::complex<double> (0, 1) *
      amp[1025] + std::complex<double> (0, 1) * amp[1023] -
      std::complex<double> (0, 1) * amp[1031] + std::complex<double> (0, 1) *
      amp[1029];
  jamp[6] = +amp[109] + amp[115] + amp[153] + amp[154] + amp[156] + amp[159] +
      amp[160] + amp[162] - std::complex<double> (0, 1) * amp[164] -
      std::complex<double> (0, 1) * amp[165] - std::complex<double> (0, 1) *
      amp[166] - std::complex<double> (0, 1) * amp[167] - amp[169] -
      std::complex<double> (0, 1) * amp[170] - amp[172] - amp[175] -
      std::complex<double> (0, 1) * amp[176] - amp[178] - std::complex<double>
      (0, 1) * amp[192] - amp[193] - std::complex<double> (0, 1) * amp[194] -
      amp[196] - std::complex<double> (0, 1) * amp[197] - std::complex<double>
      (0, 1) * amp[198] - amp[199] - std::complex<double> (0, 1) * amp[200] -
      amp[202] - std::complex<double> (0, 1) * amp[203] + amp[206] - amp[204] +
      amp[209] - amp[207] + amp[212] - amp[210] + amp[215] - amp[213] +
      amp[617] + amp[623] + std::complex<double> (0, 1) * amp[625] + amp[626] +
      std::complex<double> (0, 1) * amp[627] + amp[628] + std::complex<double>
      (0, 1) * amp[629] + std::complex<double> (0, 1) * amp[631] + amp[632] +
      std::complex<double> (0, 1) * amp[633] + amp[634] + std::complex<double>
      (0, 1) * amp[635] - std::complex<double> (0, 1) * amp[636] + amp[638] +
      amp[640] - std::complex<double> (0, 1) * amp[642] + amp[644] + amp[646] -
      amp[662] - amp[661] - amp[665] - amp[664] - amp[668] - amp[667] -
      amp[671] - amp[670] - std::complex<double> (0, 1) * amp[674] +
      std::complex<double> (0, 1) * amp[680] - std::complex<double> (0, 1) *
      amp[678] + std::complex<double> (0, 1) * amp[682] - std::complex<double>
      (0, 1) * amp[683] + amp[684] + std::complex<double> (0, 1) * amp[689] -
      std::complex<double> (0, 1) * amp[687] - std::complex<double> (0, 1) *
      amp[692] + std::complex<double> (0, 1) * amp[698] - std::complex<double>
      (0, 1) * amp[696] + std::complex<double> (0, 1) * amp[700] -
      std::complex<double> (0, 1) * amp[701] + amp[702] + std::complex<double>
      (0, 1) * amp[707] - std::complex<double> (0, 1) * amp[705] + amp[709] -
      std::complex<double> (0, 1) * amp[710] + amp[711] - std::complex<double>
      (0, 1) * amp[713] + amp[714] + amp[717] - std::complex<double> (0, 1) *
      amp[718] + amp[719] - std::complex<double> (0, 1) * amp[721] + amp[722] -
      std::complex<double> (0, 1) * amp[726] - std::complex<double> (0, 1) *
      amp[728] - amp[729] + std::complex<double> (0, 1) * amp[732] +
      std::complex<double> (0, 1) * amp[731] - std::complex<double> (0, 1) *
      amp[735] + amp[736] - std::complex<double> (0, 1) * amp[737] +
      std::complex<double> (0, 1) * amp[741] + std::complex<double> (0, 1) *
      amp[740] - std::complex<double> (0, 1) * amp[744] - std::complex<double>
      (0, 1) * amp[746] - amp[747] + std::complex<double> (0, 1) * amp[750] +
      std::complex<double> (0, 1) * amp[749] - std::complex<double> (0, 1) *
      amp[753] + amp[754] - std::complex<double> (0, 1) * amp[755] +
      std::complex<double> (0, 1) * amp[759] + std::complex<double> (0, 1) *
      amp[758] - std::complex<double> (0, 1) * amp[901] - std::complex<double>
      (0, 1) * amp[902] - amp[903] + std::complex<double> (0, 1) * amp[908] +
      std::complex<double> (0, 1) * amp[907] - std::complex<double> (0, 1) *
      amp[911] + std::complex<double> (0, 1) * amp[917] + std::complex<double>
      (0, 1) * amp[916] - std::complex<double> (0, 1) * amp[919] -
      std::complex<double> (0, 1) * amp[920] - amp[921] + std::complex<double>
      (0, 1) * amp[926] + std::complex<double> (0, 1) * amp[925] -
      std::complex<double> (0, 1) * amp[929] + std::complex<double> (0, 1) *
      amp[935] + std::complex<double> (0, 1) * amp[934] + std::complex<double>
      (0, 1) * amp[941] + std::complex<double> (0, 1) * amp[940] - amp[944] -
      amp[943] + std::complex<double> (0, 1) * amp[947] + std::complex<double>
      (0, 1) * amp[946] + std::complex<double> (0, 1) * amp[953] +
      std::complex<double> (0, 1) * amp[952] - amp[956] - amp[955] +
      std::complex<double> (0, 1) * amp[959] + std::complex<double> (0, 1) *
      amp[958] + amp[986] - amp[984] + std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[987] + std::complex<double> (0, 1) *
      amp[995] - std::complex<double> (0, 1) * amp[993] + amp[998] - amp[996] +
      std::complex<double> (0, 1) * amp[1001] - std::complex<double> (0, 1) *
      amp[999] + std::complex<double> (0, 1) * amp[1007] - std::complex<double>
      (0, 1) * amp[1005];
  jamp[7] = +amp[108] + amp[114] + amp[137] + amp[138] + amp[140] + amp[143] +
      amp[144] + amp[146] - std::complex<double> (0, 1) * amp[148] -
      std::complex<double> (0, 1) * amp[149] - std::complex<double> (0, 1) *
      amp[150] - std::complex<double> (0, 1) * amp[151] - amp[181] -
      std::complex<double> (0, 1) * amp[182] - amp[184] - amp[187] -
      std::complex<double> (0, 1) * amp[188] - amp[190] + std::complex<double>
      (0, 1) * amp[192] + amp[193] + std::complex<double> (0, 1) * amp[194] +
      amp[196] + std::complex<double> (0, 1) * amp[197] + std::complex<double>
      (0, 1) * amp[198] + amp[199] + std::complex<double> (0, 1) * amp[200] +
      amp[202] + std::complex<double> (0, 1) * amp[203] + amp[204] + amp[205] +
      amp[207] + amp[208] + amp[210] + amp[211] + amp[213] + amp[214] +
      amp[557] + amp[563] + std::complex<double> (0, 1) * amp[565] + amp[566] +
      std::complex<double> (0, 1) * amp[567] + amp[568] + std::complex<double>
      (0, 1) * amp[569] + std::complex<double> (0, 1) * amp[571] + amp[572] +
      std::complex<double> (0, 1) * amp[573] + amp[574] + std::complex<double>
      (0, 1) * amp[575] - std::complex<double> (0, 1) * amp[576] + amp[578] +
      amp[580] - std::complex<double> (0, 1) * amp[582] + amp[584] + amp[586] -
      amp[602] - amp[601] - amp[605] - amp[604] - amp[608] - amp[607] -
      amp[611] - amp[610] - std::complex<double> (0, 1) * amp[676] +
      std::complex<double> (0, 1) * amp[678] + std::complex<double> (0, 1) *
      amp[679] + std::complex<double> (0, 1) * amp[681] - std::complex<double>
      (0, 1) * amp[685] + amp[686] + std::complex<double> (0, 1) * amp[687] +
      std::complex<double> (0, 1) * amp[688] - std::complex<double> (0, 1) *
      amp[694] + std::complex<double> (0, 1) * amp[696] + std::complex<double>
      (0, 1) * amp[697] + std::complex<double> (0, 1) * amp[699] -
      std::complex<double> (0, 1) * amp[703] + amp[704] + std::complex<double>
      (0, 1) * amp[705] + std::complex<double> (0, 1) * amp[706] - amp[709] +
      std::complex<double> (0, 1) * amp[710] - amp[711] + std::complex<double>
      (0, 1) * amp[713] - amp[714] - amp[717] + std::complex<double> (0, 1) *
      amp[718] - amp[719] + std::complex<double> (0, 1) * amp[721] - amp[722] -
      std::complex<double> (0, 1) * amp[778] - std::complex<double> (0, 1) *
      amp[780] - amp[781] + std::complex<double> (0, 1) * amp[784] +
      std::complex<double> (0, 1) * amp[783] - std::complex<double> (0, 1) *
      amp[787] + amp[788] - std::complex<double> (0, 1) * amp[789] +
      std::complex<double> (0, 1) * amp[793] + std::complex<double> (0, 1) *
      amp[792] - std::complex<double> (0, 1) * amp[796] - std::complex<double>
      (0, 1) * amp[798] - amp[799] + std::complex<double> (0, 1) * amp[802] +
      std::complex<double> (0, 1) * amp[801] - std::complex<double> (0, 1) *
      amp[805] + amp[806] - std::complex<double> (0, 1) * amp[807] +
      std::complex<double> (0, 1) * amp[811] + std::complex<double> (0, 1) *
      amp[810] + std::complex<double> (0, 1) * amp[901] + std::complex<double>
      (0, 1) * amp[902] + amp[903] - std::complex<double> (0, 1) * amp[908] -
      std::complex<double> (0, 1) * amp[907] + std::complex<double> (0, 1) *
      amp[911] - std::complex<double> (0, 1) * amp[917] - std::complex<double>
      (0, 1) * amp[916] + std::complex<double> (0, 1) * amp[919] +
      std::complex<double> (0, 1) * amp[920] + amp[921] - std::complex<double>
      (0, 1) * amp[926] - std::complex<double> (0, 1) * amp[925] +
      std::complex<double> (0, 1) * amp[929] - std::complex<double> (0, 1) *
      amp[935] - std::complex<double> (0, 1) * amp[934] + std::complex<double>
      (0, 1) * amp[965] + std::complex<double> (0, 1) * amp[964] - amp[968] -
      amp[967] + std::complex<double> (0, 1) * amp[971] + std::complex<double>
      (0, 1) * amp[970] + std::complex<double> (0, 1) * amp[977] +
      std::complex<double> (0, 1) * amp[976] - amp[980] - amp[979] +
      std::complex<double> (0, 1) * amp[983] + std::complex<double> (0, 1) *
      amp[982] + amp[984] + amp[985] + std::complex<double> (0, 1) * amp[987] +
      std::complex<double> (0, 1) * amp[988] + std::complex<double> (0, 1) *
      amp[993] + std::complex<double> (0, 1) * amp[994] + amp[996] + amp[997] +
      std::complex<double> (0, 1) * amp[999] + std::complex<double> (0, 1) *
      amp[1000] + std::complex<double> (0, 1) * amp[1005] +
      std::complex<double> (0, 1) * amp[1006];
  jamp[8] = +amp[111] + amp[117] + amp[152] + amp[155] + amp[157] + amp[158] +
      amp[161] + amp[163] + std::complex<double> (0, 1) * amp[164] +
      std::complex<double> (0, 1) * amp[165] + std::complex<double> (0, 1) *
      amp[166] + std::complex<double> (0, 1) * amp[167] + amp[169] +
      std::complex<double> (0, 1) * amp[170] + amp[172] + amp[175] +
      std::complex<double> (0, 1) * amp[176] + amp[178] - std::complex<double>
      (0, 1) * amp[180] + amp[181] - std::complex<double> (0, 1) * amp[183] +
      amp[184] - std::complex<double> (0, 1) * amp[185] - std::complex<double>
      (0, 1) * amp[186] + amp[187] - std::complex<double> (0, 1) * amp[189] +
      amp[190] - std::complex<double> (0, 1) * amp[191] - amp[206] - amp[205] -
      amp[209] - amp[208] - amp[212] - amp[211] - amp[215] - amp[214] +
      amp[612] + amp[618] + std::complex<double> (0, 1) * amp[636] - amp[638] -
      amp[640] + std::complex<double> (0, 1) * amp[642] - amp[644] - amp[646] -
      std::complex<double> (0, 1) * amp[648] + amp[650] - std::complex<double>
      (0, 1) * amp[651] + amp[652] - std::complex<double> (0, 1) * amp[653] -
      std::complex<double> (0, 1) * amp[654] + amp[656] - std::complex<double>
      (0, 1) * amp[657] + amp[658] - std::complex<double> (0, 1) * amp[659] +
      amp[660] + amp[661] + amp[663] + amp[664] + amp[666] + amp[667] +
      amp[669] + amp[670] + std::complex<double> (0, 1) * amp[726] +
      std::complex<double> (0, 1) * amp[728] + amp[729] - std::complex<double>
      (0, 1) * amp[732] - std::complex<double> (0, 1) * amp[731] +
      std::complex<double> (0, 1) * amp[735] - amp[736] + std::complex<double>
      (0, 1) * amp[737] - std::complex<double> (0, 1) * amp[741] -
      std::complex<double> (0, 1) * amp[740] + std::complex<double> (0, 1) *
      amp[744] + std::complex<double> (0, 1) * amp[746] + amp[747] -
      std::complex<double> (0, 1) * amp[750] - std::complex<double> (0, 1) *
      amp[749] + std::complex<double> (0, 1) * amp[753] - amp[754] +
      std::complex<double> (0, 1) * amp[755] - std::complex<double> (0, 1) *
      amp[759] - std::complex<double> (0, 1) * amp[758] - std::complex<double>
      (0, 1) * amp[776] + std::complex<double> (0, 1) * amp[780] + amp[781] -
      std::complex<double> (0, 1) * amp[782] - std::complex<double> (0, 1) *
      amp[783] + std::complex<double> (0, 1) * amp[789] - std::complex<double>
      (0, 1) * amp[791] - std::complex<double> (0, 1) * amp[792] -
      std::complex<double> (0, 1) * amp[794] + std::complex<double> (0, 1) *
      amp[798] + amp[799] - std::complex<double> (0, 1) * amp[800] -
      std::complex<double> (0, 1) * amp[801] + std::complex<double> (0, 1) *
      amp[807] - std::complex<double> (0, 1) * amp[809] - std::complex<double>
      (0, 1) * amp[810] - amp[812] - std::complex<double> (0, 1) * amp[814] -
      amp[816] - std::complex<double> (0, 1) * amp[817] - amp[819] - amp[820] -
      std::complex<double> (0, 1) * amp[822] - amp[824] - std::complex<double>
      (0, 1) * amp[825] - amp[827] - std::complex<double> (0, 1) * amp[830] +
      std::complex<double> (0, 1) * amp[836] + std::complex<double> (0, 1) *
      amp[835] - std::complex<double> (0, 1) * amp[838] - std::complex<double>
      (0, 1) * amp[839] + amp[840] + std::complex<double> (0, 1) * amp[845] +
      std::complex<double> (0, 1) * amp[844] - std::complex<double> (0, 1) *
      amp[848] + std::complex<double> (0, 1) * amp[854] + std::complex<double>
      (0, 1) * amp[853] - std::complex<double> (0, 1) * amp[856] -
      std::complex<double> (0, 1) * amp[857] + amp[858] + std::complex<double>
      (0, 1) * amp[863] + std::complex<double> (0, 1) * amp[862] -
      std::complex<double> (0, 1) * amp[939] - std::complex<double> (0, 1) *
      amp[940] + amp[942] + amp[943] - std::complex<double> (0, 1) * amp[945] -
      std::complex<double> (0, 1) * amp[946] - std::complex<double> (0, 1) *
      amp[951] - std::complex<double> (0, 1) * amp[952] + amp[954] + amp[955] -
      std::complex<double> (0, 1) * amp[957] - std::complex<double> (0, 1) *
      amp[958] - amp[986] - amp[985] - std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[988] - std::complex<double> (0, 1) *
      amp[995] - std::complex<double> (0, 1) * amp[994] - amp[998] - amp[997] -
      std::complex<double> (0, 1) * amp[1001] - std::complex<double> (0, 1) *
      amp[1000] - std::complex<double> (0, 1) * amp[1007] -
      std::complex<double> (0, 1) * amp[1006];
  jamp[9] = +amp[110] + amp[116] + amp[121] + amp[122] + amp[124] + amp[127] +
      amp[128] + amp[130] - std::complex<double> (0, 1) * amp[132] -
      std::complex<double> (0, 1) * amp[133] - std::complex<double> (0, 1) *
      amp[134] - std::complex<double> (0, 1) * amp[135] + std::complex<double>
      (0, 1) * amp[180] - amp[181] + std::complex<double> (0, 1) * amp[183] -
      amp[184] + std::complex<double> (0, 1) * amp[185] + std::complex<double>
      (0, 1) * amp[186] - amp[187] + std::complex<double> (0, 1) * amp[189] -
      amp[190] + std::complex<double> (0, 1) * amp[191] + amp[193] -
      std::complex<double> (0, 1) * amp[195] + amp[196] + amp[199] -
      std::complex<double> (0, 1) * amp[201] + amp[202] + amp[204] + amp[205] +
      amp[207] + amp[208] + amp[210] + amp[211] + amp[213] + amp[214] +
      amp[436] + amp[442] - std::complex<double> (0, 1) * amp[445] - amp[446] -
      std::complex<double> (0, 1) * amp[447] - amp[448] - std::complex<double>
      (0, 1) * amp[449] - std::complex<double> (0, 1) * amp[451] - amp[452] -
      std::complex<double> (0, 1) * amp[453] - amp[454] - std::complex<double>
      (0, 1) * amp[455] - std::complex<double> (0, 1) * amp[469] - amp[470] -
      amp[472] - std::complex<double> (0, 1) * amp[475] - amp[476] - amp[478] +
      amp[482] - amp[480] + amp[485] - amp[483] + amp[488] - amp[486] +
      amp[491] - amp[489] + std::complex<double> (0, 1) * amp[776] -
      std::complex<double> (0, 1) * amp[780] - amp[781] + std::complex<double>
      (0, 1) * amp[782] + std::complex<double> (0, 1) * amp[783] -
      std::complex<double> (0, 1) * amp[789] + std::complex<double> (0, 1) *
      amp[791] + std::complex<double> (0, 1) * amp[792] + std::complex<double>
      (0, 1) * amp[794] - std::complex<double> (0, 1) * amp[798] - amp[799] +
      std::complex<double> (0, 1) * amp[800] + std::complex<double> (0, 1) *
      amp[801] - std::complex<double> (0, 1) * amp[807] + std::complex<double>
      (0, 1) * amp[809] + std::complex<double> (0, 1) * amp[810] + amp[812] +
      std::complex<double> (0, 1) * amp[814] + amp[816] + std::complex<double>
      (0, 1) * amp[817] + amp[819] + amp[820] + std::complex<double> (0, 1) *
      amp[822] + amp[824] + std::complex<double> (0, 1) * amp[825] + amp[827] -
      std::complex<double> (0, 1) * amp[832] - std::complex<double> (0, 1) *
      amp[834] - std::complex<double> (0, 1) * amp[835] - std::complex<double>
      (0, 1) * amp[837] - std::complex<double> (0, 1) * amp[841] - amp[842] -
      std::complex<double> (0, 1) * amp[843] - std::complex<double> (0, 1) *
      amp[844] - std::complex<double> (0, 1) * amp[850] - std::complex<double>
      (0, 1) * amp[852] - std::complex<double> (0, 1) * amp[853] -
      std::complex<double> (0, 1) * amp[855] - std::complex<double> (0, 1) *
      amp[859] - amp[860] - std::complex<double> (0, 1) * amp[861] -
      std::complex<double> (0, 1) * amp[862] + std::complex<double> (0, 1) *
      amp[902] + amp[903] - std::complex<double> (0, 1) * amp[904] -
      std::complex<double> (0, 1) * amp[906] - std::complex<double> (0, 1) *
      amp[907] + std::complex<double> (0, 1) * amp[911] - std::complex<double>
      (0, 1) * amp[913] - amp[914] - std::complex<double> (0, 1) * amp[915] -
      std::complex<double> (0, 1) * amp[916] + std::complex<double> (0, 1) *
      amp[920] + amp[921] - std::complex<double> (0, 1) * amp[922] -
      std::complex<double> (0, 1) * amp[924] - std::complex<double> (0, 1) *
      amp[925] + std::complex<double> (0, 1) * amp[929] - std::complex<double>
      (0, 1) * amp[931] - amp[932] - std::complex<double> (0, 1) * amp[933] -
      std::complex<double> (0, 1) * amp[934] + amp[984] + amp[985] +
      std::complex<double> (0, 1) * amp[987] + std::complex<double> (0, 1) *
      amp[988] + std::complex<double> (0, 1) * amp[993] + std::complex<double>
      (0, 1) * amp[994] + amp[996] + amp[997] + std::complex<double> (0, 1) *
      amp[999] + std::complex<double> (0, 1) * amp[1000] + std::complex<double>
      (0, 1) * amp[1005] + std::complex<double> (0, 1) * amp[1006] +
      std::complex<double> (0, 1) * amp[1013] - std::complex<double> (0, 1) *
      amp[1011] + amp[1016] - amp[1014] + std::complex<double> (0, 1) *
      amp[1019] - std::complex<double> (0, 1) * amp[1017] +
      std::complex<double> (0, 1) * amp[1025] - std::complex<double> (0, 1) *
      amp[1023] + amp[1028] - amp[1026] + std::complex<double> (0, 1) *
      amp[1031] - std::complex<double> (0, 1) * amp[1029];
  jamp[10] = +amp[113] + amp[119] + amp[136] + amp[139] + amp[141] + amp[142] +
      amp[145] + amp[147] + std::complex<double> (0, 1) * amp[148] +
      std::complex<double> (0, 1) * amp[149] + std::complex<double> (0, 1) *
      amp[150] + std::complex<double> (0, 1) * amp[151] - std::complex<double>
      (0, 1) * amp[168] + amp[169] - std::complex<double> (0, 1) * amp[171] +
      amp[172] - std::complex<double> (0, 1) * amp[173] - std::complex<double>
      (0, 1) * amp[174] + amp[175] - std::complex<double> (0, 1) * amp[177] +
      amp[178] - std::complex<double> (0, 1) * amp[179] + amp[181] +
      std::complex<double> (0, 1) * amp[182] + amp[184] + amp[187] +
      std::complex<double> (0, 1) * amp[188] + amp[190] - amp[206] - amp[205] -
      amp[209] - amp[208] - amp[212] - amp[211] - amp[215] - amp[214] +
      amp[552] + amp[558] + std::complex<double> (0, 1) * amp[576] - amp[578] -
      amp[580] + std::complex<double> (0, 1) * amp[582] - amp[584] - amp[586] -
      std::complex<double> (0, 1) * amp[588] + amp[590] - std::complex<double>
      (0, 1) * amp[591] + amp[592] - std::complex<double> (0, 1) * amp[593] -
      std::complex<double> (0, 1) * amp[594] + amp[596] - std::complex<double>
      (0, 1) * amp[597] + amp[598] - std::complex<double> (0, 1) * amp[599] +
      amp[600] + amp[601] + amp[603] + amp[604] + amp[606] + amp[607] +
      amp[609] + amp[610] - std::complex<double> (0, 1) * amp[724] +
      std::complex<double> (0, 1) * amp[728] + amp[729] - std::complex<double>
      (0, 1) * amp[730] - std::complex<double> (0, 1) * amp[731] +
      std::complex<double> (0, 1) * amp[737] - std::complex<double> (0, 1) *
      amp[739] - std::complex<double> (0, 1) * amp[740] - std::complex<double>
      (0, 1) * amp[742] + std::complex<double> (0, 1) * amp[746] + amp[747] -
      std::complex<double> (0, 1) * amp[748] - std::complex<double> (0, 1) *
      amp[749] + std::complex<double> (0, 1) * amp[755] - std::complex<double>
      (0, 1) * amp[757] - std::complex<double> (0, 1) * amp[758] - amp[760] -
      std::complex<double> (0, 1) * amp[762] - amp[764] - std::complex<double>
      (0, 1) * amp[765] - amp[767] - amp[768] - std::complex<double> (0, 1) *
      amp[770] - amp[772] - std::complex<double> (0, 1) * amp[773] - amp[775] +
      std::complex<double> (0, 1) * amp[778] + std::complex<double> (0, 1) *
      amp[780] + amp[781] - std::complex<double> (0, 1) * amp[784] -
      std::complex<double> (0, 1) * amp[783] + std::complex<double> (0, 1) *
      amp[787] - amp[788] + std::complex<double> (0, 1) * amp[789] -
      std::complex<double> (0, 1) * amp[793] - std::complex<double> (0, 1) *
      amp[792] + std::complex<double> (0, 1) * amp[796] + std::complex<double>
      (0, 1) * amp[798] + amp[799] - std::complex<double> (0, 1) * amp[802] -
      std::complex<double> (0, 1) * amp[801] + std::complex<double> (0, 1) *
      amp[805] - amp[806] + std::complex<double> (0, 1) * amp[807] -
      std::complex<double> (0, 1) * amp[811] - std::complex<double> (0, 1) *
      amp[810] - std::complex<double> (0, 1) * amp[866] + std::complex<double>
      (0, 1) * amp[872] + std::complex<double> (0, 1) * amp[871] -
      std::complex<double> (0, 1) * amp[874] - std::complex<double> (0, 1) *
      amp[875] + amp[876] + std::complex<double> (0, 1) * amp[881] +
      std::complex<double> (0, 1) * amp[880] - std::complex<double> (0, 1) *
      amp[884] + std::complex<double> (0, 1) * amp[890] + std::complex<double>
      (0, 1) * amp[889] - std::complex<double> (0, 1) * amp[892] -
      std::complex<double> (0, 1) * amp[893] + amp[894] + std::complex<double>
      (0, 1) * amp[899] + std::complex<double> (0, 1) * amp[898] -
      std::complex<double> (0, 1) * amp[963] - std::complex<double> (0, 1) *
      amp[964] + amp[966] + amp[967] - std::complex<double> (0, 1) * amp[969] -
      std::complex<double> (0, 1) * amp[970] - std::complex<double> (0, 1) *
      amp[975] - std::complex<double> (0, 1) * amp[976] + amp[978] + amp[979] -
      std::complex<double> (0, 1) * amp[981] - std::complex<double> (0, 1) *
      amp[982] - amp[986] - amp[985] - std::complex<double> (0, 1) * amp[989] -
      std::complex<double> (0, 1) * amp[988] - std::complex<double> (0, 1) *
      amp[995] - std::complex<double> (0, 1) * amp[994] - amp[998] - amp[997] -
      std::complex<double> (0, 1) * amp[1001] - std::complex<double> (0, 1) *
      amp[1000] - std::complex<double> (0, 1) * amp[1007] -
      std::complex<double> (0, 1) * amp[1006];
  jamp[11] = +amp[112] + amp[118] + amp[120] + amp[123] + amp[125] + amp[126] +
      amp[129] + amp[131] + std::complex<double> (0, 1) * amp[132] +
      std::complex<double> (0, 1) * amp[133] + std::complex<double> (0, 1) *
      amp[134] + std::complex<double> (0, 1) * amp[135] + std::complex<double>
      (0, 1) * amp[168] - amp[169] + std::complex<double> (0, 1) * amp[171] -
      amp[172] + std::complex<double> (0, 1) * amp[173] + std::complex<double>
      (0, 1) * amp[174] - amp[175] + std::complex<double> (0, 1) * amp[177] -
      amp[178] + std::complex<double> (0, 1) * amp[179] - amp[193] +
      std::complex<double> (0, 1) * amp[195] - amp[196] - amp[199] +
      std::complex<double> (0, 1) * amp[201] - amp[202] + amp[206] - amp[204] +
      amp[209] - amp[207] + amp[212] - amp[210] + amp[215] - amp[213] +
      amp[434] + amp[440] - std::complex<double> (0, 1) * amp[457] - amp[458] -
      std::complex<double> (0, 1) * amp[459] - amp[460] - std::complex<double>
      (0, 1) * amp[461] - std::complex<double> (0, 1) * amp[463] - amp[464] -
      std::complex<double> (0, 1) * amp[465] - amp[466] - std::complex<double>
      (0, 1) * amp[467] + std::complex<double> (0, 1) * amp[469] + amp[470] +
      amp[472] + std::complex<double> (0, 1) * amp[475] + amp[476] + amp[478] +
      amp[480] + amp[481] + amp[483] + amp[484] + amp[486] + amp[487] +
      amp[489] + amp[490] + std::complex<double> (0, 1) * amp[724] -
      std::complex<double> (0, 1) * amp[728] - amp[729] + std::complex<double>
      (0, 1) * amp[730] + std::complex<double> (0, 1) * amp[731] -
      std::complex<double> (0, 1) * amp[737] + std::complex<double> (0, 1) *
      amp[739] + std::complex<double> (0, 1) * amp[740] + std::complex<double>
      (0, 1) * amp[742] - std::complex<double> (0, 1) * amp[746] - amp[747] +
      std::complex<double> (0, 1) * amp[748] + std::complex<double> (0, 1) *
      amp[749] - std::complex<double> (0, 1) * amp[755] + std::complex<double>
      (0, 1) * amp[757] + std::complex<double> (0, 1) * amp[758] + amp[760] +
      std::complex<double> (0, 1) * amp[762] + amp[764] + std::complex<double>
      (0, 1) * amp[765] + amp[767] + amp[768] + std::complex<double> (0, 1) *
      amp[770] + amp[772] + std::complex<double> (0, 1) * amp[773] + amp[775] -
      std::complex<double> (0, 1) * amp[868] - std::complex<double> (0, 1) *
      amp[870] - std::complex<double> (0, 1) * amp[871] - std::complex<double>
      (0, 1) * amp[873] - std::complex<double> (0, 1) * amp[877] - amp[878] -
      std::complex<double> (0, 1) * amp[879] - std::complex<double> (0, 1) *
      amp[880] - std::complex<double> (0, 1) * amp[886] - std::complex<double>
      (0, 1) * amp[888] - std::complex<double> (0, 1) * amp[889] -
      std::complex<double> (0, 1) * amp[891] - std::complex<double> (0, 1) *
      amp[895] - amp[896] - std::complex<double> (0, 1) * amp[897] -
      std::complex<double> (0, 1) * amp[898] - std::complex<double> (0, 1) *
      amp[902] - amp[903] + std::complex<double> (0, 1) * amp[904] +
      std::complex<double> (0, 1) * amp[906] + std::complex<double> (0, 1) *
      amp[907] - std::complex<double> (0, 1) * amp[911] + std::complex<double>
      (0, 1) * amp[913] + amp[914] + std::complex<double> (0, 1) * amp[915] +
      std::complex<double> (0, 1) * amp[916] - std::complex<double> (0, 1) *
      amp[920] - amp[921] + std::complex<double> (0, 1) * amp[922] +
      std::complex<double> (0, 1) * amp[924] + std::complex<double> (0, 1) *
      amp[925] - std::complex<double> (0, 1) * amp[929] + std::complex<double>
      (0, 1) * amp[931] + amp[932] + std::complex<double> (0, 1) * amp[933] +
      std::complex<double> (0, 1) * amp[934] + amp[986] - amp[984] +
      std::complex<double> (0, 1) * amp[989] - std::complex<double> (0, 1) *
      amp[987] + std::complex<double> (0, 1) * amp[995] - std::complex<double>
      (0, 1) * amp[993] + amp[998] - amp[996] + std::complex<double> (0, 1) *
      amp[1001] - std::complex<double> (0, 1) * amp[999] + std::complex<double>
      (0, 1) * amp[1007] - std::complex<double> (0, 1) * amp[1005] +
      std::complex<double> (0, 1) * amp[1011] + std::complex<double> (0, 1) *
      amp[1012] + amp[1014] + amp[1015] + std::complex<double> (0, 1) *
      amp[1017] + std::complex<double> (0, 1) * amp[1018] +
      std::complex<double> (0, 1) * amp[1023] + std::complex<double> (0, 1) *
      amp[1024] + amp[1026] + amp[1027] + std::complex<double> (0, 1) *
      amp[1029] + std::complex<double> (0, 1) * amp[1030];
  jamp[12] = +amp[217] + amp[223] + amp[261] + amp[262] + amp[264] + amp[267] +
      amp[268] + amp[270] - std::complex<double> (0, 1) * amp[272] -
      std::complex<double> (0, 1) * amp[273] - std::complex<double> (0, 1) *
      amp[274] - std::complex<double> (0, 1) * amp[275] - amp[277] -
      std::complex<double> (0, 1) * amp[278] - amp[280] - amp[283] -
      std::complex<double> (0, 1) * amp[284] - amp[286] - std::complex<double>
      (0, 1) * amp[300] - amp[301] - std::complex<double> (0, 1) * amp[302] -
      amp[304] - std::complex<double> (0, 1) * amp[305] - std::complex<double>
      (0, 1) * amp[306] - amp[307] - std::complex<double> (0, 1) * amp[308] -
      amp[310] - std::complex<double> (0, 1) * amp[311] + amp[314] - amp[312] +
      amp[317] - amp[315] + amp[320] - amp[318] + amp[323] - amp[321] +
      amp[615] + amp[621] - std::complex<double> (0, 1) * amp[624] + amp[626] +
      amp[628] - std::complex<double> (0, 1) * amp[630] + amp[632] + amp[634] +
      std::complex<double> (0, 1) * amp[637] + amp[638] + std::complex<double>
      (0, 1) * amp[639] + amp[640] + std::complex<double> (0, 1) * amp[641] +
      std::complex<double> (0, 1) * amp[643] + amp[644] + std::complex<double>
      (0, 1) * amp[645] + amp[646] + std::complex<double> (0, 1) * amp[647] -
      amp[662] - amp[661] - amp[665] - amp[664] - amp[668] - amp[667] -
      amp[671] - amp[670] - std::complex<double> (0, 1) * amp[674] -
      std::complex<double> (0, 1) * amp[676] - amp[677] + std::complex<double>
      (0, 1) * amp[680] + std::complex<double> (0, 1) * amp[679] -
      std::complex<double> (0, 1) * amp[683] + amp[684] - std::complex<double>
      (0, 1) * amp[685] + std::complex<double> (0, 1) * amp[689] +
      std::complex<double> (0, 1) * amp[688] - std::complex<double> (0, 1) *
      amp[692] - std::complex<double> (0, 1) * amp[694] - amp[695] +
      std::complex<double> (0, 1) * amp[698] + std::complex<double> (0, 1) *
      amp[697] - std::complex<double> (0, 1) * amp[701] + amp[702] -
      std::complex<double> (0, 1) * amp[703] + std::complex<double> (0, 1) *
      amp[707] + std::complex<double> (0, 1) * amp[706] - std::complex<double>
      (0, 1) * amp[726] + std::complex<double> (0, 1) * amp[732] -
      std::complex<double> (0, 1) * amp[730] + std::complex<double> (0, 1) *
      amp[734] - std::complex<double> (0, 1) * amp[735] + amp[736] +
      std::complex<double> (0, 1) * amp[741] - std::complex<double> (0, 1) *
      amp[739] - std::complex<double> (0, 1) * amp[744] + std::complex<double>
      (0, 1) * amp[750] - std::complex<double> (0, 1) * amp[748] +
      std::complex<double> (0, 1) * amp[752] - std::complex<double> (0, 1) *
      amp[753] + amp[754] + std::complex<double> (0, 1) * amp[759] -
      std::complex<double> (0, 1) * amp[757] + amp[761] - std::complex<double>
      (0, 1) * amp[762] + amp[763] - std::complex<double> (0, 1) * amp[765] +
      amp[766] + amp[769] - std::complex<double> (0, 1) * amp[770] + amp[771] -
      std::complex<double> (0, 1) * amp[773] + amp[774] - std::complex<double>
      (0, 1) * amp[865] - std::complex<double> (0, 1) * amp[866] - amp[867] +
      std::complex<double> (0, 1) * amp[872] + std::complex<double> (0, 1) *
      amp[871] - std::complex<double> (0, 1) * amp[875] + std::complex<double>
      (0, 1) * amp[881] + std::complex<double> (0, 1) * amp[880] -
      std::complex<double> (0, 1) * amp[883] - std::complex<double> (0, 1) *
      amp[884] - amp[885] + std::complex<double> (0, 1) * amp[890] +
      std::complex<double> (0, 1) * amp[889] - std::complex<double> (0, 1) *
      amp[893] + std::complex<double> (0, 1) * amp[899] + std::complex<double>
      (0, 1) * amp[898] + std::complex<double> (0, 1) * amp[941] +
      std::complex<double> (0, 1) * amp[940] - amp[944] - amp[943] +
      std::complex<double> (0, 1) * amp[947] + std::complex<double> (0, 1) *
      amp[946] + std::complex<double> (0, 1) * amp[953] + std::complex<double>
      (0, 1) * amp[952] - amp[956] - amp[955] + std::complex<double> (0, 1) *
      amp[959] + std::complex<double> (0, 1) * amp[958] + amp[962] - amp[960] +
      std::complex<double> (0, 1) * amp[965] - std::complex<double> (0, 1) *
      amp[963] + std::complex<double> (0, 1) * amp[971] - std::complex<double>
      (0, 1) * amp[969] + amp[974] - amp[972] + std::complex<double> (0, 1) *
      amp[977] - std::complex<double> (0, 1) * amp[975] + std::complex<double>
      (0, 1) * amp[983] - std::complex<double> (0, 1) * amp[981];
  jamp[13] = +amp[216] + amp[222] + amp[245] + amp[246] + amp[248] + amp[251] +
      amp[252] + amp[254] - std::complex<double> (0, 1) * amp[256] -
      std::complex<double> (0, 1) * amp[257] - std::complex<double> (0, 1) *
      amp[258] - std::complex<double> (0, 1) * amp[259] - amp[289] -
      std::complex<double> (0, 1) * amp[290] - amp[292] - amp[295] -
      std::complex<double> (0, 1) * amp[296] - amp[298] + std::complex<double>
      (0, 1) * amp[300] + amp[301] + std::complex<double> (0, 1) * amp[302] +
      amp[304] + std::complex<double> (0, 1) * amp[305] + std::complex<double>
      (0, 1) * amp[306] + amp[307] + std::complex<double> (0, 1) * amp[308] +
      amp[310] + std::complex<double> (0, 1) * amp[311] + amp[312] + amp[313] +
      amp[315] + amp[316] + amp[318] + amp[319] + amp[321] + amp[322] +
      amp[497] + amp[503] + std::complex<double> (0, 1) * amp[505] + amp[506] +
      std::complex<double> (0, 1) * amp[507] + amp[508] + std::complex<double>
      (0, 1) * amp[509] + std::complex<double> (0, 1) * amp[511] + amp[512] +
      std::complex<double> (0, 1) * amp[513] + amp[514] + std::complex<double>
      (0, 1) * amp[515] - std::complex<double> (0, 1) * amp[516] + amp[518] +
      amp[520] - std::complex<double> (0, 1) * amp[522] + amp[524] + amp[526] -
      amp[542] - amp[541] - amp[545] - amp[544] - amp[548] - amp[547] -
      amp[551] - amp[550] - std::complex<double> (0, 1) * amp[728] +
      std::complex<double> (0, 1) * amp[730] + std::complex<double> (0, 1) *
      amp[731] + std::complex<double> (0, 1) * amp[733] - std::complex<double>
      (0, 1) * amp[737] + amp[738] + std::complex<double> (0, 1) * amp[739] +
      std::complex<double> (0, 1) * amp[740] - std::complex<double> (0, 1) *
      amp[746] + std::complex<double> (0, 1) * amp[748] + std::complex<double>
      (0, 1) * amp[749] + std::complex<double> (0, 1) * amp[751] -
      std::complex<double> (0, 1) * amp[755] + amp[756] + std::complex<double>
      (0, 1) * amp[757] + std::complex<double> (0, 1) * amp[758] - amp[761] +
      std::complex<double> (0, 1) * amp[762] - amp[763] + std::complex<double>
      (0, 1) * amp[765] - amp[766] - amp[769] + std::complex<double> (0, 1) *
      amp[770] - amp[771] + std::complex<double> (0, 1) * amp[773] - amp[774] -
      std::complex<double> (0, 1) * amp[778] - amp[779] - std::complex<double>
      (0, 1) * amp[780] + std::complex<double> (0, 1) * amp[784] +
      std::complex<double> (0, 1) * amp[783] - std::complex<double> (0, 1) *
      amp[787] - std::complex<double> (0, 1) * amp[789] + amp[790] +
      std::complex<double> (0, 1) * amp[793] + std::complex<double> (0, 1) *
      amp[792] - std::complex<double> (0, 1) * amp[796] - amp[797] -
      std::complex<double> (0, 1) * amp[798] + std::complex<double> (0, 1) *
      amp[802] + std::complex<double> (0, 1) * amp[801] - std::complex<double>
      (0, 1) * amp[805] - std::complex<double> (0, 1) * amp[807] + amp[808] +
      std::complex<double> (0, 1) * amp[811] + std::complex<double> (0, 1) *
      amp[810] + std::complex<double> (0, 1) * amp[865] + std::complex<double>
      (0, 1) * amp[866] + amp[867] - std::complex<double> (0, 1) * amp[872] -
      std::complex<double> (0, 1) * amp[871] + std::complex<double> (0, 1) *
      amp[875] - std::complex<double> (0, 1) * amp[881] - std::complex<double>
      (0, 1) * amp[880] + std::complex<double> (0, 1) * amp[883] +
      std::complex<double> (0, 1) * amp[884] + amp[885] - std::complex<double>
      (0, 1) * amp[890] - std::complex<double> (0, 1) * amp[889] +
      std::complex<double> (0, 1) * amp[893] - std::complex<double> (0, 1) *
      amp[899] - std::complex<double> (0, 1) * amp[898] + amp[960] + amp[961] +
      std::complex<double> (0, 1) * amp[963] + std::complex<double> (0, 1) *
      amp[964] + std::complex<double> (0, 1) * amp[969] + std::complex<double>
      (0, 1) * amp[970] + amp[972] + amp[973] + std::complex<double> (0, 1) *
      amp[975] + std::complex<double> (0, 1) * amp[976] + std::complex<double>
      (0, 1) * amp[981] + std::complex<double> (0, 1) * amp[982] +
      std::complex<double> (0, 1) * amp[989] + std::complex<double> (0, 1) *
      amp[988] - amp[992] - amp[991] + std::complex<double> (0, 1) * amp[995] +
      std::complex<double> (0, 1) * amp[994] + std::complex<double> (0, 1) *
      amp[1001] + std::complex<double> (0, 1) * amp[1000] - amp[1004] -
      amp[1003] + std::complex<double> (0, 1) * amp[1007] +
      std::complex<double> (0, 1) * amp[1006];
  jamp[14] = +amp[219] + amp[225] + amp[260] + amp[263] + amp[265] + amp[266] +
      amp[269] + amp[271] + std::complex<double> (0, 1) * amp[272] +
      std::complex<double> (0, 1) * amp[273] + std::complex<double> (0, 1) *
      amp[274] + std::complex<double> (0, 1) * amp[275] + amp[277] +
      std::complex<double> (0, 1) * amp[278] + amp[280] + amp[283] +
      std::complex<double> (0, 1) * amp[284] + amp[286] - std::complex<double>
      (0, 1) * amp[288] + amp[289] - std::complex<double> (0, 1) * amp[291] +
      amp[292] - std::complex<double> (0, 1) * amp[293] - std::complex<double>
      (0, 1) * amp[294] + amp[295] - std::complex<double> (0, 1) * amp[297] +
      amp[298] - std::complex<double> (0, 1) * amp[299] - amp[314] - amp[313] -
      amp[317] - amp[316] - amp[320] - amp[319] - amp[323] - amp[322] +
      amp[613] + amp[619] + std::complex<double> (0, 1) * amp[624] - amp[626] -
      amp[628] + std::complex<double> (0, 1) * amp[630] - amp[632] - amp[634] +
      std::complex<double> (0, 1) * amp[648] - amp[650] + std::complex<double>
      (0, 1) * amp[651] - amp[652] + std::complex<double> (0, 1) * amp[653] +
      std::complex<double> (0, 1) * amp[654] - amp[656] + std::complex<double>
      (0, 1) * amp[657] - amp[658] + std::complex<double> (0, 1) * amp[659] +
      amp[662] - amp[660] + amp[665] - amp[663] + amp[668] - amp[666] +
      amp[671] - amp[669] + std::complex<double> (0, 1) * amp[674] +
      std::complex<double> (0, 1) * amp[676] + amp[677] - std::complex<double>
      (0, 1) * amp[680] - std::complex<double> (0, 1) * amp[679] +
      std::complex<double> (0, 1) * amp[683] - amp[684] + std::complex<double>
      (0, 1) * amp[685] - std::complex<double> (0, 1) * amp[689] -
      std::complex<double> (0, 1) * amp[688] + std::complex<double> (0, 1) *
      amp[692] + std::complex<double> (0, 1) * amp[694] + amp[695] -
      std::complex<double> (0, 1) * amp[698] - std::complex<double> (0, 1) *
      amp[697] + std::complex<double> (0, 1) * amp[701] - amp[702] +
      std::complex<double> (0, 1) * amp[703] - std::complex<double> (0, 1) *
      amp[707] - std::complex<double> (0, 1) * amp[706] - std::complex<double>
      (0, 1) * amp[777] + std::complex<double> (0, 1) * amp[778] + amp[779] -
      std::complex<double> (0, 1) * amp[784] + std::complex<double> (0, 1) *
      amp[782] + std::complex<double> (0, 1) * amp[787] - std::complex<double>
      (0, 1) * amp[793] + std::complex<double> (0, 1) * amp[791] -
      std::complex<double> (0, 1) * amp[795] + std::complex<double> (0, 1) *
      amp[796] + amp[797] - std::complex<double> (0, 1) * amp[802] +
      std::complex<double> (0, 1) * amp[800] + std::complex<double> (0, 1) *
      amp[805] - std::complex<double> (0, 1) * amp[811] + std::complex<double>
      (0, 1) * amp[809] + amp[812] + std::complex<double> (0, 1) * amp[814] +
      amp[816] + std::complex<double> (0, 1) * amp[817] + amp[819] + amp[820] +
      std::complex<double> (0, 1) * amp[822] + amp[824] + std::complex<double>
      (0, 1) * amp[825] + amp[827] + std::complex<double> (0, 1) * amp[830] -
      std::complex<double> (0, 1) * amp[836] - std::complex<double> (0, 1) *
      amp[835] + std::complex<double> (0, 1) * amp[838] + std::complex<double>
      (0, 1) * amp[839] - amp[840] - std::complex<double> (0, 1) * amp[845] -
      std::complex<double> (0, 1) * amp[844] + std::complex<double> (0, 1) *
      amp[848] - std::complex<double> (0, 1) * amp[854] - std::complex<double>
      (0, 1) * amp[853] + std::complex<double> (0, 1) * amp[856] +
      std::complex<double> (0, 1) * amp[857] - amp[858] - std::complex<double>
      (0, 1) * amp[863] - std::complex<double> (0, 1) * amp[862] -
      std::complex<double> (0, 1) * amp[941] + std::complex<double> (0, 1) *
      amp[939] + amp[944] - amp[942] - std::complex<double> (0, 1) * amp[947] +
      std::complex<double> (0, 1) * amp[945] - std::complex<double> (0, 1) *
      amp[953] + std::complex<double> (0, 1) * amp[951] + amp[956] - amp[954] -
      std::complex<double> (0, 1) * amp[959] + std::complex<double> (0, 1) *
      amp[957] - amp[962] - amp[961] - std::complex<double> (0, 1) * amp[965] -
      std::complex<double> (0, 1) * amp[964] - std::complex<double> (0, 1) *
      amp[971] - std::complex<double> (0, 1) * amp[970] - amp[974] - amp[973] -
      std::complex<double> (0, 1) * amp[977] - std::complex<double> (0, 1) *
      amp[976] - std::complex<double> (0, 1) * amp[983] - std::complex<double>
      (0, 1) * amp[982];
  jamp[15] = +amp[218] + amp[224] + amp[229] + amp[230] + amp[232] + amp[235] +
      amp[236] + amp[238] - std::complex<double> (0, 1) * amp[240] -
      std::complex<double> (0, 1) * amp[241] - std::complex<double> (0, 1) *
      amp[242] - std::complex<double> (0, 1) * amp[243] + std::complex<double>
      (0, 1) * amp[288] - amp[289] + std::complex<double> (0, 1) * amp[291] -
      amp[292] + std::complex<double> (0, 1) * amp[293] + std::complex<double>
      (0, 1) * amp[294] - amp[295] + std::complex<double> (0, 1) * amp[297] -
      amp[298] + std::complex<double> (0, 1) * amp[299] + amp[301] -
      std::complex<double> (0, 1) * amp[303] + amp[304] + amp[307] -
      std::complex<double> (0, 1) * amp[309] + amp[310] + amp[312] + amp[313] +
      amp[315] + amp[316] + amp[318] + amp[319] + amp[321] + amp[322] +
      amp[437] + amp[443] + std::complex<double> (0, 1) * amp[445] + amp[446] +
      std::complex<double> (0, 1) * amp[447] + amp[448] + std::complex<double>
      (0, 1) * amp[449] + std::complex<double> (0, 1) * amp[451] + amp[452] +
      std::complex<double> (0, 1) * amp[453] + amp[454] + std::complex<double>
      (0, 1) * amp[455] - std::complex<double> (0, 1) * amp[456] + amp[458] +
      amp[460] - std::complex<double> (0, 1) * amp[462] + amp[464] + amp[466] -
      amp[482] - amp[481] - amp[485] - amp[484] - amp[488] - amp[487] -
      amp[491] - amp[490] + std::complex<double> (0, 1) * amp[777] -
      std::complex<double> (0, 1) * amp[778] - amp[779] + std::complex<double>
      (0, 1) * amp[784] - std::complex<double> (0, 1) * amp[782] -
      std::complex<double> (0, 1) * amp[787] + std::complex<double> (0, 1) *
      amp[793] - std::complex<double> (0, 1) * amp[791] + std::complex<double>
      (0, 1) * amp[795] - std::complex<double> (0, 1) * amp[796] - amp[797] +
      std::complex<double> (0, 1) * amp[802] - std::complex<double> (0, 1) *
      amp[800] - std::complex<double> (0, 1) * amp[805] + std::complex<double>
      (0, 1) * amp[811] - std::complex<double> (0, 1) * amp[809] - amp[812] -
      std::complex<double> (0, 1) * amp[814] - amp[816] - std::complex<double>
      (0, 1) * amp[817] - amp[819] - amp[820] - std::complex<double> (0, 1) *
      amp[822] - amp[824] - std::complex<double> (0, 1) * amp[825] - amp[827] +
      std::complex<double> (0, 1) * amp[832] + std::complex<double> (0, 1) *
      amp[834] + std::complex<double> (0, 1) * amp[835] + std::complex<double>
      (0, 1) * amp[837] + std::complex<double> (0, 1) * amp[841] + amp[842] +
      std::complex<double> (0, 1) * amp[843] + std::complex<double> (0, 1) *
      amp[844] + std::complex<double> (0, 1) * amp[850] + std::complex<double>
      (0, 1) * amp[852] + std::complex<double> (0, 1) * amp[853] +
      std::complex<double> (0, 1) * amp[855] + std::complex<double> (0, 1) *
      amp[859] + amp[860] + std::complex<double> (0, 1) * amp[861] +
      std::complex<double> (0, 1) * amp[862] + std::complex<double> (0, 1) *
      amp[866] + amp[867] + std::complex<double> (0, 1) * amp[868] -
      std::complex<double> (0, 1) * amp[872] + std::complex<double> (0, 1) *
      amp[870] + std::complex<double> (0, 1) * amp[875] + std::complex<double>
      (0, 1) * amp[877] + amp[878] - std::complex<double> (0, 1) * amp[881] +
      std::complex<double> (0, 1) * amp[879] + std::complex<double> (0, 1) *
      amp[884] + amp[885] + std::complex<double> (0, 1) * amp[886] -
      std::complex<double> (0, 1) * amp[890] + std::complex<double> (0, 1) *
      amp[888] + std::complex<double> (0, 1) * amp[893] + std::complex<double>
      (0, 1) * amp[895] + amp[896] - std::complex<double> (0, 1) * amp[899] +
      std::complex<double> (0, 1) * amp[897] + amp[960] + amp[961] +
      std::complex<double> (0, 1) * amp[963] + std::complex<double> (0, 1) *
      amp[964] + std::complex<double> (0, 1) * amp[969] + std::complex<double>
      (0, 1) * amp[970] + amp[972] + amp[973] + std::complex<double> (0, 1) *
      amp[975] + std::complex<double> (0, 1) * amp[976] + std::complex<double>
      (0, 1) * amp[981] + std::complex<double> (0, 1) * amp[982] -
      std::complex<double> (0, 1) * amp[1013] - std::complex<double> (0, 1) *
      amp[1012] - amp[1016] - amp[1015] - std::complex<double> (0, 1) *
      amp[1019] - std::complex<double> (0, 1) * amp[1018] -
      std::complex<double> (0, 1) * amp[1025] - std::complex<double> (0, 1) *
      amp[1024] - amp[1028] - amp[1027] - std::complex<double> (0, 1) *
      amp[1031] - std::complex<double> (0, 1) * amp[1030];
  jamp[16] = +amp[221] + amp[227] + amp[244] + amp[247] + amp[249] + amp[250] +
      amp[253] + amp[255] + std::complex<double> (0, 1) * amp[256] +
      std::complex<double> (0, 1) * amp[257] + std::complex<double> (0, 1) *
      amp[258] + std::complex<double> (0, 1) * amp[259] - std::complex<double>
      (0, 1) * amp[276] + amp[277] - std::complex<double> (0, 1) * amp[279] +
      amp[280] - std::complex<double> (0, 1) * amp[281] - std::complex<double>
      (0, 1) * amp[282] + amp[283] - std::complex<double> (0, 1) * amp[285] +
      amp[286] - std::complex<double> (0, 1) * amp[287] + amp[289] +
      std::complex<double> (0, 1) * amp[290] + amp[292] + amp[295] +
      std::complex<double> (0, 1) * amp[296] + amp[298] - amp[314] - amp[313] -
      amp[317] - amp[316] - amp[320] - amp[319] - amp[323] - amp[322] +
      amp[492] + amp[498] + std::complex<double> (0, 1) * amp[516] - amp[518] -
      amp[520] + std::complex<double> (0, 1) * amp[522] - amp[524] - amp[526] -
      std::complex<double> (0, 1) * amp[528] + amp[530] - std::complex<double>
      (0, 1) * amp[531] + amp[532] - std::complex<double> (0, 1) * amp[533] -
      std::complex<double> (0, 1) * amp[534] + amp[536] - std::complex<double>
      (0, 1) * amp[537] + amp[538] - std::complex<double> (0, 1) * amp[539] +
      amp[540] + amp[541] + amp[543] + amp[544] + amp[546] + amp[547] +
      amp[549] + amp[550] - std::complex<double> (0, 1) * amp[672] +
      std::complex<double> (0, 1) * amp[676] + amp[677] - std::complex<double>
      (0, 1) * amp[678] - std::complex<double> (0, 1) * amp[679] +
      std::complex<double> (0, 1) * amp[685] - std::complex<double> (0, 1) *
      amp[687] - std::complex<double> (0, 1) * amp[688] - std::complex<double>
      (0, 1) * amp[690] + std::complex<double> (0, 1) * amp[694] + amp[695] -
      std::complex<double> (0, 1) * amp[696] - std::complex<double> (0, 1) *
      amp[697] + std::complex<double> (0, 1) * amp[703] - std::complex<double>
      (0, 1) * amp[705] - std::complex<double> (0, 1) * amp[706] - amp[708] -
      std::complex<double> (0, 1) * amp[710] - amp[712] - std::complex<double>
      (0, 1) * amp[713] - amp[715] - amp[716] - std::complex<double> (0, 1) *
      amp[718] - amp[720] - std::complex<double> (0, 1) * amp[721] - amp[723] +
      std::complex<double> (0, 1) * amp[778] + amp[779] + std::complex<double>
      (0, 1) * amp[780] - std::complex<double> (0, 1) * amp[784] -
      std::complex<double> (0, 1) * amp[783] + std::complex<double> (0, 1) *
      amp[787] + std::complex<double> (0, 1) * amp[789] - amp[790] -
      std::complex<double> (0, 1) * amp[793] - std::complex<double> (0, 1) *
      amp[792] + std::complex<double> (0, 1) * amp[796] + amp[797] +
      std::complex<double> (0, 1) * amp[798] - std::complex<double> (0, 1) *
      amp[802] - std::complex<double> (0, 1) * amp[801] + std::complex<double>
      (0, 1) * amp[805] + std::complex<double> (0, 1) * amp[807] - amp[808] -
      std::complex<double> (0, 1) * amp[811] - std::complex<double> (0, 1) *
      amp[810] - std::complex<double> (0, 1) * amp[902] + std::complex<double>
      (0, 1) * amp[908] + std::complex<double> (0, 1) * amp[907] -
      std::complex<double> (0, 1) * amp[910] - std::complex<double> (0, 1) *
      amp[911] + amp[912] + std::complex<double> (0, 1) * amp[917] +
      std::complex<double> (0, 1) * amp[916] - std::complex<double> (0, 1) *
      amp[920] + std::complex<double> (0, 1) * amp[926] + std::complex<double>
      (0, 1) * amp[925] - std::complex<double> (0, 1) * amp[928] -
      std::complex<double> (0, 1) * amp[929] + amp[930] + std::complex<double>
      (0, 1) * amp[935] + std::complex<double> (0, 1) * amp[934] - amp[962] -
      amp[961] - std::complex<double> (0, 1) * amp[965] - std::complex<double>
      (0, 1) * amp[964] - std::complex<double> (0, 1) * amp[971] -
      std::complex<double> (0, 1) * amp[970] - amp[974] - amp[973] -
      std::complex<double> (0, 1) * amp[977] - std::complex<double> (0, 1) *
      amp[976] - std::complex<double> (0, 1) * amp[983] - std::complex<double>
      (0, 1) * amp[982] - std::complex<double> (0, 1) * amp[987] -
      std::complex<double> (0, 1) * amp[988] + amp[990] + amp[991] -
      std::complex<double> (0, 1) * amp[993] - std::complex<double> (0, 1) *
      amp[994] - std::complex<double> (0, 1) * amp[999] - std::complex<double>
      (0, 1) * amp[1000] + amp[1002] + amp[1003] - std::complex<double> (0, 1)
      * amp[1005] - std::complex<double> (0, 1) * amp[1006];
  jamp[17] = +amp[220] + amp[226] + amp[228] + amp[231] + amp[233] + amp[234] +
      amp[237] + amp[239] + std::complex<double> (0, 1) * amp[240] +
      std::complex<double> (0, 1) * amp[241] + std::complex<double> (0, 1) *
      amp[242] + std::complex<double> (0, 1) * amp[243] + std::complex<double>
      (0, 1) * amp[276] - amp[277] + std::complex<double> (0, 1) * amp[279] -
      amp[280] + std::complex<double> (0, 1) * amp[281] + std::complex<double>
      (0, 1) * amp[282] - amp[283] + std::complex<double> (0, 1) * amp[285] -
      amp[286] + std::complex<double> (0, 1) * amp[287] - amp[301] +
      std::complex<double> (0, 1) * amp[303] - amp[304] - amp[307] +
      std::complex<double> (0, 1) * amp[309] - amp[310] + amp[314] - amp[312] +
      amp[317] - amp[315] + amp[320] - amp[318] + amp[323] - amp[321] +
      amp[432] + amp[438] + std::complex<double> (0, 1) * amp[456] - amp[458] -
      amp[460] + std::complex<double> (0, 1) * amp[462] - amp[464] - amp[466] -
      std::complex<double> (0, 1) * amp[468] + amp[470] - std::complex<double>
      (0, 1) * amp[471] + amp[472] - std::complex<double> (0, 1) * amp[473] -
      std::complex<double> (0, 1) * amp[474] + amp[476] - std::complex<double>
      (0, 1) * amp[477] + amp[478] - std::complex<double> (0, 1) * amp[479] +
      amp[480] + amp[481] + amp[483] + amp[484] + amp[486] + amp[487] +
      amp[489] + amp[490] + std::complex<double> (0, 1) * amp[672] -
      std::complex<double> (0, 1) * amp[676] - amp[677] + std::complex<double>
      (0, 1) * amp[678] + std::complex<double> (0, 1) * amp[679] -
      std::complex<double> (0, 1) * amp[685] + std::complex<double> (0, 1) *
      amp[687] + std::complex<double> (0, 1) * amp[688] + std::complex<double>
      (0, 1) * amp[690] - std::complex<double> (0, 1) * amp[694] - amp[695] +
      std::complex<double> (0, 1) * amp[696] + std::complex<double> (0, 1) *
      amp[697] - std::complex<double> (0, 1) * amp[703] + std::complex<double>
      (0, 1) * amp[705] + std::complex<double> (0, 1) * amp[706] + amp[708] +
      std::complex<double> (0, 1) * amp[710] + amp[712] + std::complex<double>
      (0, 1) * amp[713] + amp[715] + amp[716] + std::complex<double> (0, 1) *
      amp[718] + amp[720] + std::complex<double> (0, 1) * amp[721] + amp[723] -
      std::complex<double> (0, 1) * amp[866] - amp[867] - std::complex<double>
      (0, 1) * amp[868] + std::complex<double> (0, 1) * amp[872] -
      std::complex<double> (0, 1) * amp[870] - std::complex<double> (0, 1) *
      amp[875] - std::complex<double> (0, 1) * amp[877] - amp[878] +
      std::complex<double> (0, 1) * amp[881] - std::complex<double> (0, 1) *
      amp[879] - std::complex<double> (0, 1) * amp[884] - amp[885] -
      std::complex<double> (0, 1) * amp[886] + std::complex<double> (0, 1) *
      amp[890] - std::complex<double> (0, 1) * amp[888] - std::complex<double>
      (0, 1) * amp[893] - std::complex<double> (0, 1) * amp[895] - amp[896] +
      std::complex<double> (0, 1) * amp[899] - std::complex<double> (0, 1) *
      amp[897] + std::complex<double> (0, 1) * amp[904] - std::complex<double>
      (0, 1) * amp[908] + std::complex<double> (0, 1) * amp[906] -
      std::complex<double> (0, 1) * amp[909] + std::complex<double> (0, 1) *
      amp[913] + amp[914] - std::complex<double> (0, 1) * amp[917] +
      std::complex<double> (0, 1) * amp[915] + std::complex<double> (0, 1) *
      amp[922] - std::complex<double> (0, 1) * amp[926] + std::complex<double>
      (0, 1) * amp[924] - std::complex<double> (0, 1) * amp[927] +
      std::complex<double> (0, 1) * amp[931] + amp[932] - std::complex<double>
      (0, 1) * amp[935] + std::complex<double> (0, 1) * amp[933] + amp[962] -
      amp[960] + std::complex<double> (0, 1) * amp[965] - std::complex<double>
      (0, 1) * amp[963] + std::complex<double> (0, 1) * amp[971] -
      std::complex<double> (0, 1) * amp[969] + amp[974] - amp[972] +
      std::complex<double> (0, 1) * amp[977] - std::complex<double> (0, 1) *
      amp[975] + std::complex<double> (0, 1) * amp[983] - std::complex<double>
      (0, 1) * amp[981] + std::complex<double> (0, 1) * amp[1011] +
      std::complex<double> (0, 1) * amp[1012] + amp[1014] + amp[1015] +
      std::complex<double> (0, 1) * amp[1017] + std::complex<double> (0, 1) *
      amp[1018] + std::complex<double> (0, 1) * amp[1023] +
      std::complex<double> (0, 1) * amp[1024] + amp[1026] + amp[1027] +
      std::complex<double> (0, 1) * amp[1029] + std::complex<double> (0, 1) *
      amp[1030];
  jamp[18] = +amp[325] + amp[331] + amp[369] + amp[370] + amp[372] + amp[375] +
      amp[376] + amp[378] - std::complex<double> (0, 1) * amp[380] -
      std::complex<double> (0, 1) * amp[381] - std::complex<double> (0, 1) *
      amp[382] - std::complex<double> (0, 1) * amp[383] - amp[385] -
      std::complex<double> (0, 1) * amp[386] - amp[388] - amp[391] -
      std::complex<double> (0, 1) * amp[392] - amp[394] - std::complex<double>
      (0, 1) * amp[408] - amp[409] - std::complex<double> (0, 1) * amp[410] -
      amp[412] - std::complex<double> (0, 1) * amp[413] - std::complex<double>
      (0, 1) * amp[414] - amp[415] - std::complex<double> (0, 1) * amp[416] -
      amp[418] - std::complex<double> (0, 1) * amp[419] + amp[422] - amp[420] +
      amp[425] - amp[423] + amp[428] - amp[426] + amp[431] - amp[429] +
      amp[555] + amp[561] - std::complex<double> (0, 1) * amp[564] + amp[566] +
      amp[568] - std::complex<double> (0, 1) * amp[570] + amp[572] + amp[574] +
      std::complex<double> (0, 1) * amp[577] + amp[578] + std::complex<double>
      (0, 1) * amp[579] + amp[580] + std::complex<double> (0, 1) * amp[581] +
      std::complex<double> (0, 1) * amp[583] + amp[584] + std::complex<double>
      (0, 1) * amp[585] + amp[586] + std::complex<double> (0, 1) * amp[587] -
      amp[602] - amp[601] - amp[605] - amp[604] - amp[608] - amp[607] -
      amp[611] - amp[610] - std::complex<double> (0, 1) * amp[674] - amp[675] -
      std::complex<double> (0, 1) * amp[676] + std::complex<double> (0, 1) *
      amp[680] + std::complex<double> (0, 1) * amp[679] - std::complex<double>
      (0, 1) * amp[683] - std::complex<double> (0, 1) * amp[685] + amp[686] +
      std::complex<double> (0, 1) * amp[689] + std::complex<double> (0, 1) *
      amp[688] - std::complex<double> (0, 1) * amp[692] - amp[693] -
      std::complex<double> (0, 1) * amp[694] + std::complex<double> (0, 1) *
      amp[698] + std::complex<double> (0, 1) * amp[697] - std::complex<double>
      (0, 1) * amp[701] - std::complex<double> (0, 1) * amp[703] + amp[704] +
      std::complex<double> (0, 1) * amp[707] + std::complex<double> (0, 1) *
      amp[706] - std::complex<double> (0, 1) * amp[778] + std::complex<double>
      (0, 1) * amp[784] - std::complex<double> (0, 1) * amp[782] +
      std::complex<double> (0, 1) * amp[786] - std::complex<double> (0, 1) *
      amp[787] + amp[788] + std::complex<double> (0, 1) * amp[793] -
      std::complex<double> (0, 1) * amp[791] - std::complex<double> (0, 1) *
      amp[796] + std::complex<double> (0, 1) * amp[802] - std::complex<double>
      (0, 1) * amp[800] + std::complex<double> (0, 1) * amp[804] -
      std::complex<double> (0, 1) * amp[805] + amp[806] + std::complex<double>
      (0, 1) * amp[811] - std::complex<double> (0, 1) * amp[809] + amp[813] -
      std::complex<double> (0, 1) * amp[814] + amp[815] - std::complex<double>
      (0, 1) * amp[817] + amp[818] + amp[821] - std::complex<double> (0, 1) *
      amp[822] + amp[823] - std::complex<double> (0, 1) * amp[825] + amp[826] -
      std::complex<double> (0, 1) * amp[829] - std::complex<double> (0, 1) *
      amp[830] - amp[831] + std::complex<double> (0, 1) * amp[836] +
      std::complex<double> (0, 1) * amp[835] - std::complex<double> (0, 1) *
      amp[839] + std::complex<double> (0, 1) * amp[845] + std::complex<double>
      (0, 1) * amp[844] - std::complex<double> (0, 1) * amp[847] -
      std::complex<double> (0, 1) * amp[848] - amp[849] + std::complex<double>
      (0, 1) * amp[854] + std::complex<double> (0, 1) * amp[853] -
      std::complex<double> (0, 1) * amp[857] + std::complex<double> (0, 1) *
      amp[863] + std::complex<double> (0, 1) * amp[862] + amp[938] - amp[936] +
      std::complex<double> (0, 1) * amp[941] - std::complex<double> (0, 1) *
      amp[939] + std::complex<double> (0, 1) * amp[947] - std::complex<double>
      (0, 1) * amp[945] + amp[950] - amp[948] + std::complex<double> (0, 1) *
      amp[953] - std::complex<double> (0, 1) * amp[951] + std::complex<double>
      (0, 1) * amp[959] - std::complex<double> (0, 1) * amp[957] +
      std::complex<double> (0, 1) * amp[965] + std::complex<double> (0, 1) *
      amp[964] - amp[968] - amp[967] + std::complex<double> (0, 1) * amp[971] +
      std::complex<double> (0, 1) * amp[970] + std::complex<double> (0, 1) *
      amp[977] + std::complex<double> (0, 1) * amp[976] - amp[980] - amp[979] +
      std::complex<double> (0, 1) * amp[983] + std::complex<double> (0, 1) *
      amp[982];
  jamp[19] = +amp[324] + amp[330] + amp[353] + amp[354] + amp[356] + amp[359] +
      amp[360] + amp[362] - std::complex<double> (0, 1) * amp[364] -
      std::complex<double> (0, 1) * amp[365] - std::complex<double> (0, 1) *
      amp[366] - std::complex<double> (0, 1) * amp[367] - amp[397] -
      std::complex<double> (0, 1) * amp[398] - amp[400] - amp[403] -
      std::complex<double> (0, 1) * amp[404] - amp[406] + std::complex<double>
      (0, 1) * amp[408] + amp[409] + std::complex<double> (0, 1) * amp[410] +
      amp[412] + std::complex<double> (0, 1) * amp[413] + std::complex<double>
      (0, 1) * amp[414] + amp[415] + std::complex<double> (0, 1) * amp[416] +
      amp[418] + std::complex<double> (0, 1) * amp[419] + amp[420] + amp[421] +
      amp[423] + amp[424] + amp[426] + amp[427] + amp[429] + amp[430] +
      amp[495] + amp[501] - std::complex<double> (0, 1) * amp[504] + amp[506] +
      amp[508] - std::complex<double> (0, 1) * amp[510] + amp[512] + amp[514] +
      std::complex<double> (0, 1) * amp[517] + amp[518] + std::complex<double>
      (0, 1) * amp[519] + amp[520] + std::complex<double> (0, 1) * amp[521] +
      std::complex<double> (0, 1) * amp[523] + amp[524] + std::complex<double>
      (0, 1) * amp[525] + amp[526] + std::complex<double> (0, 1) * amp[527] -
      amp[542] - amp[541] - amp[545] - amp[544] - amp[548] - amp[547] -
      amp[551] - amp[550] - std::complex<double> (0, 1) * amp[726] - amp[727] -
      std::complex<double> (0, 1) * amp[728] + std::complex<double> (0, 1) *
      amp[732] + std::complex<double> (0, 1) * amp[731] - std::complex<double>
      (0, 1) * amp[735] - std::complex<double> (0, 1) * amp[737] + amp[738] +
      std::complex<double> (0, 1) * amp[741] + std::complex<double> (0, 1) *
      amp[740] - std::complex<double> (0, 1) * amp[744] - amp[745] -
      std::complex<double> (0, 1) * amp[746] + std::complex<double> (0, 1) *
      amp[750] + std::complex<double> (0, 1) * amp[749] - std::complex<double>
      (0, 1) * amp[753] - std::complex<double> (0, 1) * amp[755] + amp[756] +
      std::complex<double> (0, 1) * amp[759] + std::complex<double> (0, 1) *
      amp[758] - std::complex<double> (0, 1) * amp[780] + std::complex<double>
      (0, 1) * amp[782] + std::complex<double> (0, 1) * amp[783] +
      std::complex<double> (0, 1) * amp[785] - std::complex<double> (0, 1) *
      amp[789] + amp[790] + std::complex<double> (0, 1) * amp[791] +
      std::complex<double> (0, 1) * amp[792] - std::complex<double> (0, 1) *
      amp[798] + std::complex<double> (0, 1) * amp[800] + std::complex<double>
      (0, 1) * amp[801] + std::complex<double> (0, 1) * amp[803] -
      std::complex<double> (0, 1) * amp[807] + amp[808] + std::complex<double>
      (0, 1) * amp[809] + std::complex<double> (0, 1) * amp[810] - amp[813] +
      std::complex<double> (0, 1) * amp[814] - amp[815] + std::complex<double>
      (0, 1) * amp[817] - amp[818] - amp[821] + std::complex<double> (0, 1) *
      amp[822] - amp[823] + std::complex<double> (0, 1) * amp[825] - amp[826] +
      std::complex<double> (0, 1) * amp[829] + std::complex<double> (0, 1) *
      amp[830] + amp[831] - std::complex<double> (0, 1) * amp[836] -
      std::complex<double> (0, 1) * amp[835] + std::complex<double> (0, 1) *
      amp[839] - std::complex<double> (0, 1) * amp[845] - std::complex<double>
      (0, 1) * amp[844] + std::complex<double> (0, 1) * amp[847] +
      std::complex<double> (0, 1) * amp[848] + amp[849] - std::complex<double>
      (0, 1) * amp[854] - std::complex<double> (0, 1) * amp[853] +
      std::complex<double> (0, 1) * amp[857] - std::complex<double> (0, 1) *
      amp[863] - std::complex<double> (0, 1) * amp[862] + amp[936] + amp[937] +
      std::complex<double> (0, 1) * amp[939] + std::complex<double> (0, 1) *
      amp[940] + std::complex<double> (0, 1) * amp[945] + std::complex<double>
      (0, 1) * amp[946] + amp[948] + amp[949] + std::complex<double> (0, 1) *
      amp[951] + std::complex<double> (0, 1) * amp[952] + std::complex<double>
      (0, 1) * amp[957] + std::complex<double> (0, 1) * amp[958] +
      std::complex<double> (0, 1) * amp[989] + std::complex<double> (0, 1) *
      amp[988] - amp[992] - amp[991] + std::complex<double> (0, 1) * amp[995] +
      std::complex<double> (0, 1) * amp[994] + std::complex<double> (0, 1) *
      amp[1001] + std::complex<double> (0, 1) * amp[1000] - amp[1004] -
      amp[1003] + std::complex<double> (0, 1) * amp[1007] +
      std::complex<double> (0, 1) * amp[1006];
  jamp[20] = +amp[327] + amp[333] + amp[368] + amp[371] + amp[373] + amp[374] +
      amp[377] + amp[379] + std::complex<double> (0, 1) * amp[380] +
      std::complex<double> (0, 1) * amp[381] + std::complex<double> (0, 1) *
      amp[382] + std::complex<double> (0, 1) * amp[383] + amp[385] +
      std::complex<double> (0, 1) * amp[386] + amp[388] + amp[391] +
      std::complex<double> (0, 1) * amp[392] + amp[394] - std::complex<double>
      (0, 1) * amp[396] + amp[397] - std::complex<double> (0, 1) * amp[399] +
      amp[400] - std::complex<double> (0, 1) * amp[401] - std::complex<double>
      (0, 1) * amp[402] + amp[403] - std::complex<double> (0, 1) * amp[405] +
      amp[406] - std::complex<double> (0, 1) * amp[407] - amp[422] - amp[421] -
      amp[425] - amp[424] - amp[428] - amp[427] - amp[431] - amp[430] +
      amp[553] + amp[559] + std::complex<double> (0, 1) * amp[564] - amp[566] -
      amp[568] + std::complex<double> (0, 1) * amp[570] - amp[572] - amp[574] +
      std::complex<double> (0, 1) * amp[588] - amp[590] + std::complex<double>
      (0, 1) * amp[591] - amp[592] + std::complex<double> (0, 1) * amp[593] +
      std::complex<double> (0, 1) * amp[594] - amp[596] + std::complex<double>
      (0, 1) * amp[597] - amp[598] + std::complex<double> (0, 1) * amp[599] +
      amp[602] - amp[600] + amp[605] - amp[603] + amp[608] - amp[606] +
      amp[611] - amp[609] + std::complex<double> (0, 1) * amp[674] + amp[675] +
      std::complex<double> (0, 1) * amp[676] - std::complex<double> (0, 1) *
      amp[680] - std::complex<double> (0, 1) * amp[679] + std::complex<double>
      (0, 1) * amp[683] + std::complex<double> (0, 1) * amp[685] - amp[686] -
      std::complex<double> (0, 1) * amp[689] - std::complex<double> (0, 1) *
      amp[688] + std::complex<double> (0, 1) * amp[692] + amp[693] +
      std::complex<double> (0, 1) * amp[694] - std::complex<double> (0, 1) *
      amp[698] - std::complex<double> (0, 1) * amp[697] + std::complex<double>
      (0, 1) * amp[701] + std::complex<double> (0, 1) * amp[703] - amp[704] -
      std::complex<double> (0, 1) * amp[707] - std::complex<double> (0, 1) *
      amp[706] - std::complex<double> (0, 1) * amp[725] + std::complex<double>
      (0, 1) * amp[726] + amp[727] - std::complex<double> (0, 1) * amp[732] +
      std::complex<double> (0, 1) * amp[730] + std::complex<double> (0, 1) *
      amp[735] - std::complex<double> (0, 1) * amp[741] + std::complex<double>
      (0, 1) * amp[739] - std::complex<double> (0, 1) * amp[743] +
      std::complex<double> (0, 1) * amp[744] + amp[745] - std::complex<double>
      (0, 1) * amp[750] + std::complex<double> (0, 1) * amp[748] +
      std::complex<double> (0, 1) * amp[753] - std::complex<double> (0, 1) *
      amp[759] + std::complex<double> (0, 1) * amp[757] + amp[760] +
      std::complex<double> (0, 1) * amp[762] + amp[764] + std::complex<double>
      (0, 1) * amp[765] + amp[767] + amp[768] + std::complex<double> (0, 1) *
      amp[770] + amp[772] + std::complex<double> (0, 1) * amp[773] + amp[775] +
      std::complex<double> (0, 1) * amp[866] - std::complex<double> (0, 1) *
      amp[872] - std::complex<double> (0, 1) * amp[871] + std::complex<double>
      (0, 1) * amp[874] + std::complex<double> (0, 1) * amp[875] - amp[876] -
      std::complex<double> (0, 1) * amp[881] - std::complex<double> (0, 1) *
      amp[880] + std::complex<double> (0, 1) * amp[884] - std::complex<double>
      (0, 1) * amp[890] - std::complex<double> (0, 1) * amp[889] +
      std::complex<double> (0, 1) * amp[892] + std::complex<double> (0, 1) *
      amp[893] - amp[894] - std::complex<double> (0, 1) * amp[899] -
      std::complex<double> (0, 1) * amp[898] - amp[938] - amp[937] -
      std::complex<double> (0, 1) * amp[941] - std::complex<double> (0, 1) *
      amp[940] - std::complex<double> (0, 1) * amp[947] - std::complex<double>
      (0, 1) * amp[946] - amp[950] - amp[949] - std::complex<double> (0, 1) *
      amp[953] - std::complex<double> (0, 1) * amp[952] - std::complex<double>
      (0, 1) * amp[959] - std::complex<double> (0, 1) * amp[958] -
      std::complex<double> (0, 1) * amp[965] + std::complex<double> (0, 1) *
      amp[963] + amp[968] - amp[966] - std::complex<double> (0, 1) * amp[971] +
      std::complex<double> (0, 1) * amp[969] - std::complex<double> (0, 1) *
      amp[977] + std::complex<double> (0, 1) * amp[975] + amp[980] - amp[978] -
      std::complex<double> (0, 1) * amp[983] + std::complex<double> (0, 1) *
      amp[981];
  jamp[21] = +amp[326] + amp[332] + amp[337] + amp[338] + amp[340] + amp[343] +
      amp[344] + amp[346] - std::complex<double> (0, 1) * amp[348] -
      std::complex<double> (0, 1) * amp[349] - std::complex<double> (0, 1) *
      amp[350] - std::complex<double> (0, 1) * amp[351] + std::complex<double>
      (0, 1) * amp[396] - amp[397] + std::complex<double> (0, 1) * amp[399] -
      amp[400] + std::complex<double> (0, 1) * amp[401] + std::complex<double>
      (0, 1) * amp[402] - amp[403] + std::complex<double> (0, 1) * amp[405] -
      amp[406] + std::complex<double> (0, 1) * amp[407] + amp[409] -
      std::complex<double> (0, 1) * amp[411] + amp[412] + amp[415] -
      std::complex<double> (0, 1) * amp[417] + amp[418] + amp[420] + amp[421] +
      amp[423] + amp[424] + amp[426] + amp[427] + amp[429] + amp[430] +
      amp[435] + amp[441] - std::complex<double> (0, 1) * amp[444] + amp[446] +
      amp[448] - std::complex<double> (0, 1) * amp[450] + amp[452] + amp[454] +
      std::complex<double> (0, 1) * amp[457] + amp[458] + std::complex<double>
      (0, 1) * amp[459] + amp[460] + std::complex<double> (0, 1) * amp[461] +
      std::complex<double> (0, 1) * amp[463] + amp[464] + std::complex<double>
      (0, 1) * amp[465] + amp[466] + std::complex<double> (0, 1) * amp[467] -
      amp[482] - amp[481] - amp[485] - amp[484] - amp[488] - amp[487] -
      amp[491] - amp[490] + std::complex<double> (0, 1) * amp[725] -
      std::complex<double> (0, 1) * amp[726] - amp[727] + std::complex<double>
      (0, 1) * amp[732] - std::complex<double> (0, 1) * amp[730] -
      std::complex<double> (0, 1) * amp[735] + std::complex<double> (0, 1) *
      amp[741] - std::complex<double> (0, 1) * amp[739] + std::complex<double>
      (0, 1) * amp[743] - std::complex<double> (0, 1) * amp[744] - amp[745] +
      std::complex<double> (0, 1) * amp[750] - std::complex<double> (0, 1) *
      amp[748] - std::complex<double> (0, 1) * amp[753] + std::complex<double>
      (0, 1) * amp[759] - std::complex<double> (0, 1) * amp[757] - amp[760] -
      std::complex<double> (0, 1) * amp[762] - amp[764] - std::complex<double>
      (0, 1) * amp[765] - amp[767] - amp[768] - std::complex<double> (0, 1) *
      amp[770] - amp[772] - std::complex<double> (0, 1) * amp[773] - amp[775] +
      std::complex<double> (0, 1) * amp[830] + amp[831] + std::complex<double>
      (0, 1) * amp[832] - std::complex<double> (0, 1) * amp[836] +
      std::complex<double> (0, 1) * amp[834] + std::complex<double> (0, 1) *
      amp[839] + std::complex<double> (0, 1) * amp[841] + amp[842] -
      std::complex<double> (0, 1) * amp[845] + std::complex<double> (0, 1) *
      amp[843] + std::complex<double> (0, 1) * amp[848] + amp[849] +
      std::complex<double> (0, 1) * amp[850] - std::complex<double> (0, 1) *
      amp[854] + std::complex<double> (0, 1) * amp[852] + std::complex<double>
      (0, 1) * amp[857] + std::complex<double> (0, 1) * amp[859] + amp[860] -
      std::complex<double> (0, 1) * amp[863] + std::complex<double> (0, 1) *
      amp[861] + std::complex<double> (0, 1) * amp[868] + std::complex<double>
      (0, 1) * amp[870] + std::complex<double> (0, 1) * amp[871] +
      std::complex<double> (0, 1) * amp[873] + std::complex<double> (0, 1) *
      amp[877] + amp[878] + std::complex<double> (0, 1) * amp[879] +
      std::complex<double> (0, 1) * amp[880] + std::complex<double> (0, 1) *
      amp[886] + std::complex<double> (0, 1) * amp[888] + std::complex<double>
      (0, 1) * amp[889] + std::complex<double> (0, 1) * amp[891] +
      std::complex<double> (0, 1) * amp[895] + amp[896] + std::complex<double>
      (0, 1) * amp[897] + std::complex<double> (0, 1) * amp[898] + amp[936] +
      amp[937] + std::complex<double> (0, 1) * amp[939] + std::complex<double>
      (0, 1) * amp[940] + std::complex<double> (0, 1) * amp[945] +
      std::complex<double> (0, 1) * amp[946] + amp[948] + amp[949] +
      std::complex<double> (0, 1) * amp[951] + std::complex<double> (0, 1) *
      amp[952] + std::complex<double> (0, 1) * amp[957] + std::complex<double>
      (0, 1) * amp[958] - std::complex<double> (0, 1) * amp[1013] -
      std::complex<double> (0, 1) * amp[1012] - amp[1016] - amp[1015] -
      std::complex<double> (0, 1) * amp[1019] - std::complex<double> (0, 1) *
      amp[1018] - std::complex<double> (0, 1) * amp[1025] -
      std::complex<double> (0, 1) * amp[1024] - amp[1028] - amp[1027] -
      std::complex<double> (0, 1) * amp[1031] - std::complex<double> (0, 1) *
      amp[1030];
  jamp[22] = +amp[329] + amp[335] + amp[352] + amp[355] + amp[357] + amp[358] +
      amp[361] + amp[363] + std::complex<double> (0, 1) * amp[364] +
      std::complex<double> (0, 1) * amp[365] + std::complex<double> (0, 1) *
      amp[366] + std::complex<double> (0, 1) * amp[367] - std::complex<double>
      (0, 1) * amp[384] + amp[385] - std::complex<double> (0, 1) * amp[387] +
      amp[388] - std::complex<double> (0, 1) * amp[389] - std::complex<double>
      (0, 1) * amp[390] + amp[391] - std::complex<double> (0, 1) * amp[393] +
      amp[394] - std::complex<double> (0, 1) * amp[395] + amp[397] +
      std::complex<double> (0, 1) * amp[398] + amp[400] + amp[403] +
      std::complex<double> (0, 1) * amp[404] + amp[406] - amp[422] - amp[421] -
      amp[425] - amp[424] - amp[428] - amp[427] - amp[431] - amp[430] +
      amp[493] + amp[499] + std::complex<double> (0, 1) * amp[504] - amp[506] -
      amp[508] + std::complex<double> (0, 1) * amp[510] - amp[512] - amp[514] +
      std::complex<double> (0, 1) * amp[528] - amp[530] + std::complex<double>
      (0, 1) * amp[531] - amp[532] + std::complex<double> (0, 1) * amp[533] +
      std::complex<double> (0, 1) * amp[534] - amp[536] + std::complex<double>
      (0, 1) * amp[537] - amp[538] + std::complex<double> (0, 1) * amp[539] +
      amp[542] - amp[540] + amp[545] - amp[543] + amp[548] - amp[546] +
      amp[551] - amp[549] - std::complex<double> (0, 1) * amp[673] +
      std::complex<double> (0, 1) * amp[674] + amp[675] - std::complex<double>
      (0, 1) * amp[680] + std::complex<double> (0, 1) * amp[678] +
      std::complex<double> (0, 1) * amp[683] - std::complex<double> (0, 1) *
      amp[689] + std::complex<double> (0, 1) * amp[687] - std::complex<double>
      (0, 1) * amp[691] + std::complex<double> (0, 1) * amp[692] + amp[693] -
      std::complex<double> (0, 1) * amp[698] + std::complex<double> (0, 1) *
      amp[696] + std::complex<double> (0, 1) * amp[701] - std::complex<double>
      (0, 1) * amp[707] + std::complex<double> (0, 1) * amp[705] + amp[708] +
      std::complex<double> (0, 1) * amp[710] + amp[712] + std::complex<double>
      (0, 1) * amp[713] + amp[715] + amp[716] + std::complex<double> (0, 1) *
      amp[718] + amp[720] + std::complex<double> (0, 1) * amp[721] + amp[723] +
      std::complex<double> (0, 1) * amp[726] + amp[727] + std::complex<double>
      (0, 1) * amp[728] - std::complex<double> (0, 1) * amp[732] -
      std::complex<double> (0, 1) * amp[731] + std::complex<double> (0, 1) *
      amp[735] + std::complex<double> (0, 1) * amp[737] - amp[738] -
      std::complex<double> (0, 1) * amp[741] - std::complex<double> (0, 1) *
      amp[740] + std::complex<double> (0, 1) * amp[744] + amp[745] +
      std::complex<double> (0, 1) * amp[746] - std::complex<double> (0, 1) *
      amp[750] - std::complex<double> (0, 1) * amp[749] + std::complex<double>
      (0, 1) * amp[753] + std::complex<double> (0, 1) * amp[755] - amp[756] -
      std::complex<double> (0, 1) * amp[759] - std::complex<double> (0, 1) *
      amp[758] + std::complex<double> (0, 1) * amp[902] - std::complex<double>
      (0, 1) * amp[908] - std::complex<double> (0, 1) * amp[907] +
      std::complex<double> (0, 1) * amp[910] + std::complex<double> (0, 1) *
      amp[911] - amp[912] - std::complex<double> (0, 1) * amp[917] -
      std::complex<double> (0, 1) * amp[916] + std::complex<double> (0, 1) *
      amp[920] - std::complex<double> (0, 1) * amp[926] - std::complex<double>
      (0, 1) * amp[925] + std::complex<double> (0, 1) * amp[928] +
      std::complex<double> (0, 1) * amp[929] - amp[930] - std::complex<double>
      (0, 1) * amp[935] - std::complex<double> (0, 1) * amp[934] - amp[938] -
      amp[937] - std::complex<double> (0, 1) * amp[941] - std::complex<double>
      (0, 1) * amp[940] - std::complex<double> (0, 1) * amp[947] -
      std::complex<double> (0, 1) * amp[946] - amp[950] - amp[949] -
      std::complex<double> (0, 1) * amp[953] - std::complex<double> (0, 1) *
      amp[952] - std::complex<double> (0, 1) * amp[959] - std::complex<double>
      (0, 1) * amp[958] - std::complex<double> (0, 1) * amp[989] +
      std::complex<double> (0, 1) * amp[987] + amp[992] - amp[990] -
      std::complex<double> (0, 1) * amp[995] + std::complex<double> (0, 1) *
      amp[993] - std::complex<double> (0, 1) * amp[1001] + std::complex<double>
      (0, 1) * amp[999] + amp[1004] - amp[1002] - std::complex<double> (0, 1) *
      amp[1007] + std::complex<double> (0, 1) * amp[1005];
  jamp[23] = +amp[328] + amp[334] + amp[336] + amp[339] + amp[341] + amp[342] +
      amp[345] + amp[347] + std::complex<double> (0, 1) * amp[348] +
      std::complex<double> (0, 1) * amp[349] + std::complex<double> (0, 1) *
      amp[350] + std::complex<double> (0, 1) * amp[351] + std::complex<double>
      (0, 1) * amp[384] - amp[385] + std::complex<double> (0, 1) * amp[387] -
      amp[388] + std::complex<double> (0, 1) * amp[389] + std::complex<double>
      (0, 1) * amp[390] - amp[391] + std::complex<double> (0, 1) * amp[393] -
      amp[394] + std::complex<double> (0, 1) * amp[395] - amp[409] +
      std::complex<double> (0, 1) * amp[411] - amp[412] - amp[415] +
      std::complex<double> (0, 1) * amp[417] - amp[418] + amp[422] - amp[420] +
      amp[425] - amp[423] + amp[428] - amp[426] + amp[431] - amp[429] +
      amp[433] + amp[439] + std::complex<double> (0, 1) * amp[444] - amp[446] -
      amp[448] + std::complex<double> (0, 1) * amp[450] - amp[452] - amp[454] +
      std::complex<double> (0, 1) * amp[468] - amp[470] + std::complex<double>
      (0, 1) * amp[471] - amp[472] + std::complex<double> (0, 1) * amp[473] +
      std::complex<double> (0, 1) * amp[474] - amp[476] + std::complex<double>
      (0, 1) * amp[477] - amp[478] + std::complex<double> (0, 1) * amp[479] +
      amp[482] - amp[480] + amp[485] - amp[483] + amp[488] - amp[486] +
      amp[491] - amp[489] + std::complex<double> (0, 1) * amp[673] -
      std::complex<double> (0, 1) * amp[674] - amp[675] + std::complex<double>
      (0, 1) * amp[680] - std::complex<double> (0, 1) * amp[678] -
      std::complex<double> (0, 1) * amp[683] + std::complex<double> (0, 1) *
      amp[689] - std::complex<double> (0, 1) * amp[687] + std::complex<double>
      (0, 1) * amp[691] - std::complex<double> (0, 1) * amp[692] - amp[693] +
      std::complex<double> (0, 1) * amp[698] - std::complex<double> (0, 1) *
      amp[696] - std::complex<double> (0, 1) * amp[701] + std::complex<double>
      (0, 1) * amp[707] - std::complex<double> (0, 1) * amp[705] - amp[708] -
      std::complex<double> (0, 1) * amp[710] - amp[712] - std::complex<double>
      (0, 1) * amp[713] - amp[715] - amp[716] - std::complex<double> (0, 1) *
      amp[718] - amp[720] - std::complex<double> (0, 1) * amp[721] - amp[723] -
      std::complex<double> (0, 1) * amp[830] - amp[831] - std::complex<double>
      (0, 1) * amp[832] + std::complex<double> (0, 1) * amp[836] -
      std::complex<double> (0, 1) * amp[834] - std::complex<double> (0, 1) *
      amp[839] - std::complex<double> (0, 1) * amp[841] - amp[842] +
      std::complex<double> (0, 1) * amp[845] - std::complex<double> (0, 1) *
      amp[843] - std::complex<double> (0, 1) * amp[848] - amp[849] -
      std::complex<double> (0, 1) * amp[850] + std::complex<double> (0, 1) *
      amp[854] - std::complex<double> (0, 1) * amp[852] - std::complex<double>
      (0, 1) * amp[857] - std::complex<double> (0, 1) * amp[859] - amp[860] +
      std::complex<double> (0, 1) * amp[863] - std::complex<double> (0, 1) *
      amp[861] - std::complex<double> (0, 1) * amp[904] + std::complex<double>
      (0, 1) * amp[908] - std::complex<double> (0, 1) * amp[906] +
      std::complex<double> (0, 1) * amp[909] - std::complex<double> (0, 1) *
      amp[913] - amp[914] + std::complex<double> (0, 1) * amp[917] -
      std::complex<double> (0, 1) * amp[915] - std::complex<double> (0, 1) *
      amp[922] + std::complex<double> (0, 1) * amp[926] - std::complex<double>
      (0, 1) * amp[924] + std::complex<double> (0, 1) * amp[927] -
      std::complex<double> (0, 1) * amp[931] - amp[932] + std::complex<double>
      (0, 1) * amp[935] - std::complex<double> (0, 1) * amp[933] + amp[938] -
      amp[936] + std::complex<double> (0, 1) * amp[941] - std::complex<double>
      (0, 1) * amp[939] + std::complex<double> (0, 1) * amp[947] -
      std::complex<double> (0, 1) * amp[945] + amp[950] - amp[948] +
      std::complex<double> (0, 1) * amp[953] - std::complex<double> (0, 1) *
      amp[951] + std::complex<double> (0, 1) * amp[959] - std::complex<double>
      (0, 1) * amp[957] + std::complex<double> (0, 1) * amp[1013] -
      std::complex<double> (0, 1) * amp[1011] + amp[1016] - amp[1014] +
      std::complex<double> (0, 1) * amp[1019] - std::complex<double> (0, 1) *
      amp[1017] + std::complex<double> (0, 1) * amp[1025] -
      std::complex<double> (0, 1) * amp[1023] + amp[1028] - amp[1026] +
      std::complex<double> (0, 1) * amp[1031] - std::complex<double> (0, 1) *
      amp[1029];

  // Store the leading color flows for choice of color
  for(int i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return -1.;
}


double eeuugggg::get_jamp2(int i)
{
  return jamp2[0][i];
}

int eeuugggg::colorstring(int i, int j)
{
  static const int res[24][8] = {
    {5, 6, 7, 8, 3, 4, 0, 0},
    {5, 6, 8, 7, 3, 4, 0, 0},
    {5, 7, 6, 8, 3, 4, 0, 0},
    {5, 7, 8, 6, 3, 4, 0, 0},
    {5, 8, 6, 7, 3, 4, 0, 0},
    {5, 8, 7, 6, 3, 4, 0, 0},
    {6, 5, 7, 8, 3, 4, 0, 0},
    {6, 5, 8, 7, 3, 4, 0, 0},
    {6, 7, 5, 8, 3, 4, 0, 0},
    {6, 7, 8, 5, 3, 4, 0, 0},
    {6, 8, 5, 7, 3, 4, 0, 0},
    {6, 8, 7, 5, 3, 4, 0, 0},
    {7, 5, 6, 8, 3, 4, 0, 0},
    {7, 5, 8, 6, 3, 4, 0, 0},
    {7, 6, 5, 8, 3, 4, 0, 0},
    {7, 6, 8, 5, 3, 4, 0, 0},
    {7, 8, 5, 6, 3, 4, 0, 0},
    {7, 8, 6, 5, 3, 4, 0, 0},
    {8, 5, 6, 7, 3, 4, 0, 0},
    {8, 5, 7, 6, 3, 4, 0, 0},
    {8, 6, 5, 7, 3, 4, 0, 0},
    {8, 6, 7, 5, 3, 4, 0, 0},
    {8, 7, 5, 6, 3, 4, 0, 0},
    {8, 7, 6, 5, 3, 4, 0, 0}};
  return res[i][j];
}


int eeuugggg::NCol()
{
  static const int ncolor = 24;
  return ncolor;
}












