//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_sm.h"

// Initialize static instance
Parameters_sm * Parameters_sm::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_sm * Parameters_sm::getInstance()
{
  if (instance == 0)
    instance = new Parameters_sm(); 

  return instance; 
}

void Parameters_sm::setIndependentCouplings()
{
  GC_1 =  -(ee * complexi)/3.; 
  GC_2 =   (2. * ee * complexi)/3.; 
  GC_3 =  -(ee * complexi); 
  GC_50 = -(cw * ee * complexi)/(2. * sw); 
  GC_51 =  (cw * ee * complexi)/(2. * sw); 
  GC_58 = -(ee * complexi * sw)/(6. * cw); 
  GC_59 =  (ee * complexi * sw)/(2. * cw); 
}

void Parameters_sm::setDependentParameters()
{
  sqrt__aS  = sqrt(aS); 
  G         = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
}

void Parameters_sm::setDependentCouplings()
{
  GC_10 = -G; 
  GC_11 = complexi * G;
  GC_12 = complexi * G__exp__2; 
}

void Parameters_sm::setIndependentParameters(map<string, double> &MGParams)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  WH = MGParams.find("WH") ->second;
  WW = MGParams.find("WW") ->second;
  WZ = MGParams.find("WZ") ->second;
  WT = MGParams.find("WT") ->second;
  Gf = MGParams.find("GF") ->second;
  MH = MGParams.find("MH") ->second;
  MZ = MGParams.find("MZ") ->second;
  MW = MGParams.find("MW") ->second;
  MTA =MGParams.find("MTA") ->second;
  MT = MGParams.find("MT") ->second;
  MB = MGParams.find("MB") ->second;
  aS = MGParams.find("aS") ->second;
  aEWM1 = MGParams.find("aEWM1") ->second;
  ymtau = MTA;
  ymt   = MT;
  ymb   = MT;
  cw    = MW/MZ;
  cw__exp__2 = pow(cw, 2.); 
  sw = sqrt(1. - cw__exp__2);
  sw__exp_2 = pow(sw, 2.);
  conjg__CKM1x1 = 1.; 
  conjg__CKM3x3 = 1.; 
  CKM3x3 = 1.; 
  complexi = std::complex<double> (0., 1.); 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  MH__exp__2 = pow(MH, 2.); 
  aEW = 1./aEWM1; 
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.);
  g1 = ee/cw; 
  gw = ee/sw; 
  vev = (2. * MW * sw)/ee; 
  vev__exp__2 = pow(vev, 2.); 
  lam = MH__exp__2/(2. * vev__exp__2); 
  yb = (ymb * sqrt__2)/vev; 
  yt = (ymt * sqrt__2)/vev; 
  ytau = (ymtau * sqrt__2)/vev; 
  muH = sqrt(lam * vev__exp__2); 
  I1x33 = yb * conjg__CKM3x3; 
  I2x33 = yt * conjg__CKM3x3; 
  I3x33 = CKM3x3 * yt; 
  I4x33 = CKM3x3 * yb; 
  ee__exp__2 = pow(ee, 2.); 
  sw__exp__2 = pow(sw, 2.); 
  cw__exp__2 = pow(cw, 2.); 
}
