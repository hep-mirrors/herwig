// -*- C++ -*-
//
// ISGW2FormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGW2FormFactor class.
//

#include "ISGW2FormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ISGW2FormFactor::ISGW2FormFactor() { 
  // include the a(omega) piece (off by default)
  _includeaW = false;
  // values of a_S at matching scale
  _alphamuQM=0.6;
  // the quark masses
  _mdown    = 0.33*GeV;
  _mup      = 0.33*GeV;
  _mstrange = 0.55*GeV;
  _mcharm   = 1.82*GeV;
  _mbottom  = 5.20*GeV;
  // 1 1S0 parameters
  // the wavefunction parameters
  _beta1S0ud = 0.41*GeV;
  _beta1S0us = 0.44*GeV;
  _beta1S0ss = 0.53*GeV;
  _beta1S0cu = 0.45*GeV;
  _beta1S0cs = 0.56*GeV;
  _beta1S0ub = 0.43*GeV;
  _beta1S0sb = 0.54*GeV;
  _beta1S0cc = 0.88*GeV;
  _beta1S0bc = 0.92*GeV;
  // 1 3S1 parameters
  // the wavefunction parameters
  _beta3S1ud = 0.30*GeV;
  _beta3S1us = 0.33*GeV;
  _beta3S1ss = 0.37*GeV;
  _beta3S1cu = 0.38*GeV;
  _beta3S1cs = 0.44*GeV;
  _beta3S1ub = 0.40*GeV;
  _beta3S1sb = 0.49*GeV;
  _beta3S1cc = 0.62*GeV;
  _beta3S1bc = 0.75*GeV;
  // 1P parameters
  // the wavefunction parameters
  _beta1Pud = 0.28*GeV;
  _beta1Pus = 0.30*GeV;
  _beta1Pss = 0.33*GeV;
  _beta1Pcu = 0.33*GeV;
  _beta1Pcs = 0.38*GeV;
  _beta1Pub = 0.35*GeV;
  _beta1Psb = 0.41*GeV;
  _beta1Pcc = 0.52*GeV;
  _beta1Pbc = 0.60*GeV;
  // relativistic correction factors
  _CfDrho     = 0.889;
  _CfDKstar   = 0.928;
  _CfDsKstar  = 0.873;
  _CfDsphi    = 0.911;
  _CfBrho     = 0.905;
  _CfBDstar   = 0.989;
  _CfBsKstar  = 0.892;
  _CfBsDstar  = 0.984;
  _CfBcDstar  = 0.868;
  _CfBcpsi    = 0.967;
  _CfBcBsstar = 1.0;
  _CfBcBstar  = 1.0;
  // eta-eta' mixing angle
  _thetaeta = -Constants::pi/9.;
  // B_c to d cbar
  addFormFactor(-541,-411  ,0,-4, 5, 1);
  addFormFactor(-541,-413  ,1,-4, 5, 1);
  addFormFactor(-541,-415  ,2,-4, 5, 1);
  addFormFactor(-541, 10413,1,-2, 5, 1);
  addFormFactor(-541,-20413,1,-4, 5, 1);
  addFormFactor(-541, 10411,0, 4, 5, 1);
  // B_c to u cbar
  addFormFactor(-541,-421  ,0,-4, 5, 2);
  addFormFactor(-541,-423  ,1,-4, 5, 2);
  addFormFactor(-541,-425  ,2,-4, 5, 2);
  addFormFactor(-541,-10423,1,-4, 5, 2);
  addFormFactor(-541,-20423,1,-4, 5, 2);
  addFormFactor(-541,-10421,0,-4, 5, 2);
  // B_c to s cbar
  addFormFactor(-541,-431  ,0,-4, 5, 3);
  addFormFactor(-541,-433  ,1,-4, 5, 3);
  addFormFactor(-541,-435  ,2,-4, 5, 3);
  addFormFactor(-541,-10433,1,-4, 5, 3);
  addFormFactor(-541,-20433,1,-4, 5, 3);
  addFormFactor(-541,-10431,0, 4, 5, 3);
  // B_c decays to c cbar
  addFormFactor(-541, 441  ,0,-4, 5, 4);
  addFormFactor(-541, 443  ,1,-4, 5, 4);
  addFormFactor(-541, 445  ,2,-4, 5, 4);
  addFormFactor(-541, 10443,1,-4, 5, 4);
  addFormFactor(-541, 20443,1,-4, 5, 4);
  addFormFactor(-541, 10441,0, 4, 5, 4);
  // B_c to b dbar
  addFormFactor( 541, 511  ,0, 5,-4,-1);
  addFormFactor( 541, 513  ,1, 5,-4,-1);
  addFormFactor( 541, 515  ,2, 5,-4,-1);
  addFormFactor( 541, 10513,1, 5,-4,-1); 
  addFormFactor( 541, 20513,1, 5,-4,-1);
  addFormFactor( 541, 10511,0, 5,-4,-1);
  // B_c to b ubar
  addFormFactor( 541, 521  ,0, 5,-4,-2);
  addFormFactor( 541, 523  ,1, 5,-4,-2);
  addFormFactor( 541, 525  ,2, 5,-4,-2);
  addFormFactor( 541, 10523,1, 5,-4,-2); 
  addFormFactor( 541, 20523,1, 5,-4,-2);
  addFormFactor( 541, 10521,0, 5,-4,-2);
  // B_c decays to s bbar
  addFormFactor( 541, 531  ,0, 5,-4,-3);
  addFormFactor( 541, 533  ,1, 5,-4,-3);
  addFormFactor( 541, 535  ,2, 5,-4,-3);
  addFormFactor( 541, 10533,1, 5,-4,-3);
  addFormFactor( 541, 20533,1, 5,-4,-3);
  addFormFactor( 541, 10531,0, 5,-4,-3);
  // B_s to d sbar
  addFormFactor(-531, 311  ,0, 3,-5,-1);
  addFormFactor(-531, 313  ,1, 3,-5,-1);
  addFormFactor(-531, 315  ,2, 3,-5,-1);
  addFormFactor(-531, 10313,1, 3,-5,-1);
  addFormFactor(-531, 20313,1, 3,-5,-1);
  addFormFactor(-531, 10311,0, 3,-5,-1);
  // B_s to u sbar
  addFormFactor(-531, 321  ,0, 3,-5,-2);
  addFormFactor(-531, 323  ,1, 3,-5,-2);
  addFormFactor(-531, 325  ,2, 3,-5,-2);
  addFormFactor(-531, 10323,1, 3,-5,-2);
  addFormFactor(-531, 20323,1, 3,-5,-2);
  addFormFactor(-531, 10321,0, 3,-5,-2);
  // B_s to s sbar
  addFormFactor(-531, 221  ,0, 3,-5,-3);
  addFormFactor(-531, 331  ,0, 3,-5,-3);
  addFormFactor(-531, 333  ,1, 3,-5,-3);
  addFormFactor(-531, 335  ,2, 3,-5,-3);
  addFormFactor(-531, 10333,1, 3,-5,-3);
  addFormFactor(-531, 20333,1, 3,-5,-3);
  addFormFactor(-531, 10331,0, 3,-5,-3);
  // B_s decays to c sbar
  addFormFactor(-531, 431  ,0, 3,-5,-4);
  addFormFactor(-531, 433  ,1, 3,-5,-4);
  addFormFactor(-531, 435  ,2, 3,-5,-4);
  addFormFactor(-531, 10433,1, 3,-5,-4);
  addFormFactor(-531, 20433,1, 3,-5,-4);
  addFormFactor(-531, 10431,0, 3,-5,-4);
  // B_u decays to d ubar
  addFormFactor(-521,-211  ,0,-2, 5, 1);
  addFormFactor(-521,-213  ,1,-2, 5, 1);
  addFormFactor(-521,-215  ,2,-2, 5, 1);
  addFormFactor(-521,-10213,1,-2, 5, 1);
  addFormFactor(-521,-20213,1,-2, 5, 1);
  addFormFactor(-521,-10211,0,-2, 5, 1);
  // B_u to uu (I=0)
  addFormFactor(-521, 221  ,0,-2, 5, 2);
  addFormFactor(-521, 331  ,0,-2, 5, 2);
  addFormFactor(-521, 223  ,1,-2, 5, 2);
  addFormFactor(-521, 225  ,2,-2, 5, 2);
  addFormFactor(-521, 10223,1,-2, 5, 2);
  addFormFactor(-521, 20223,1,-2, 5, 2);
  addFormFactor(-521, 10221,0,-2, 5, 2);
  // B_u to uu (I=1)
  addFormFactor(-521, 111  ,0,-2, 5, 2);
  addFormFactor(-521, 113  ,1,-2, 5, 2);
  addFormFactor(-521, 115  ,2,-2, 5, 2);
  addFormFactor(-521, 10113,1,-2, 5, 2);
  addFormFactor(-521, 20113,1,-2, 5, 2);
  addFormFactor(-521, 10111,0,-2, 5, 2);
  // B_u decays to s ubar
  addFormFactor(-521,-321  ,0,-2, 5, 3);
  addFormFactor(-521,-323  ,1,-2, 5, 3);
  addFormFactor(-521,-325  ,2,-2, 5, 3);
  addFormFactor(-521,-10323,1,-2, 5, 3);
  addFormFactor(-521,-20323,1,-2, 5, 3);
  addFormFactor(-521,-10321,0,-2, 5, 3);
  // B_u decays to c ubar
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-521, 425  ,2,-2, 5, 4);
  addFormFactor(-521, 10423,1,-2, 5, 4);
  addFormFactor(-521, 20423,1,-2, 5, 4);
  addFormFactor(-521, 10421,0,-2, 5, 4);
  // B_d decays to d dbar (I=0)
  addFormFactor(-511, 221  ,0, 1,-5,-1); 
  addFormFactor(-511, 331  ,0, 1,-5,-1); 
  addFormFactor(-511, 223  ,1, 1,-5,-1); 
  addFormFactor(-511, 225  ,2, 1,-5,-1); 
  addFormFactor(-511, 10223,1, 1,-5,-1); 
  addFormFactor(-511, 20223,1, 1,-5,-1);
  addFormFactor(-511, 10221,0, 1,-5,-1);
  // B_d decays to d dbar (I=1)
  addFormFactor(-511, 111  ,0, 1,-5,-1); 
  addFormFactor(-511, 113  ,1, 1,-5,-1); 
  addFormFactor(-511, 115  ,2, 1,-5,-1); 
  addFormFactor(-511, 10113,1, 1,-5,-1); 
  addFormFactor(-511, 20113,1, 1,-5,-1);
  addFormFactor(-511, 10111,0, 1,-5,-1);
  // B_d decays to u dbar
  addFormFactor(-511, 211  ,0, 1,-5,-2); 
  addFormFactor(-511, 213  ,1, 1,-5,-2); 
  addFormFactor(-511, 215  ,2, 1,-5,-2); 
  addFormFactor(-511, 10213,1, 1,-5,-2); 
  addFormFactor(-511, 20213,1, 1,-5,-2);
  addFormFactor(-511, 10211,0, 1,-5,-2);
  // B_d decays to  s dbar
  addFormFactor(-511, 311  ,0, 1,-5,-3);
  addFormFactor(-511, 313  ,1, 1,-5,-3);
  addFormFactor(-511, 315  ,2, 1,-5,-3);
  addFormFactor(-511, 10313,1, 1,-5,-3);
  addFormFactor(-511, 20313,1, 1,-5,-3);
  addFormFactor(-511, 10311,0, 1,-5,-3);
  // B_d decays to  c dbar
  addFormFactor(-511, 411  ,0, 1,-5,-4);
  addFormFactor(-511, 413  ,1, 1,-5,-4);
  addFormFactor(-511, 415  ,2, 1,-5,-4);
  addFormFactor(-511, 10413,1, 1,-5,-4);
  addFormFactor(-511, 20413,1, 1,-5,-4);
  addFormFactor(-511, 10411,0, 1,-5,-4);
  // D_s to d sbar
  addFormFactor( 431,   311,0,-3, 4, 1);
  addFormFactor( 431,   313,1,-3, 4, 1);
  addFormFactor( 431,   315,2,-3, 4, 1);
  addFormFactor( 431, 10313,1,-3, 4, 1);
  addFormFactor( 431, 20313,1,-3, 4, 1);
  addFormFactor( 431, 10311,0,-3, 4, 1);
  // D_s to u sbar
  addFormFactor( 431,   321,0,-3, 4, 2);
  addFormFactor( 431,   323,1,-3, 4, 2);
  addFormFactor( 431,   325,2,-3, 4, 2);
  addFormFactor( 431, 10323,1,-3, 4, 2);
  addFormFactor( 431, 20323,1,-3, 4, 2);
  addFormFactor( 431, 10321,0,-3, 4, 2);
  // D_s to s sbar
  addFormFactor( 431, 221  ,0,-3, 4, 3);
  addFormFactor( 431, 331  ,0,-3, 4, 3);
  addFormFactor( 431, 333  ,1,-3, 4, 3);
  addFormFactor( 431, 335  ,2,-3, 4, 3);
  addFormFactor( 431, 10333,1,-3, 4, 3);
  addFormFactor( 431, 20333,1,-3, 4, 3);
  addFormFactor( 431, 10331,0,-3, 4, 3);
  // D0 to d ubar
  addFormFactor( 421,-211  ,0,-2, 4, 1);
  addFormFactor( 421,-213  ,1,-2, 4, 1);
  addFormFactor( 421,-215  ,2,-2, 4, 1);
  addFormFactor( 421,-10213,1,-2, 4, 1);
  addFormFactor( 421,-20213,1,-2, 4, 1);
  addFormFactor( 421,-10211,0,-2, 4, 1);
  // D0 to u ubar (I=1)
  addFormFactor( 421, 111  ,0,-2, 4, 2);
  addFormFactor( 421, 113  ,1,-2, 4, 2);
  addFormFactor( 421, 115  ,2,-2, 4, 2);
  addFormFactor( 421, 10113,1,-2, 4, 2);
  addFormFactor( 421, 20113,1,-2, 4, 2);
  addFormFactor( 421, 10111,0,-2, 4, 2);
  // D0 to u ubar (I=0)
  addFormFactor( 421, 221  ,0,-2, 4, 2);
  addFormFactor( 421, 331  ,0,-2, 4, 2);
  addFormFactor( 421, 223  ,1,-2, 4, 2);
  addFormFactor( 421, 225  ,2,-2, 4, 2);
  addFormFactor( 421, 10223,1,-2, 4, 2);
  addFormFactor( 421, 20223,1,-2, 4, 2);
  addFormFactor( 421, 10221,0,-2, 4, 2);
  // D0 to s ubar
  addFormFactor( 421,-321  ,0,-2, 4, 3);
  addFormFactor( 421,-323  ,1,-2, 4, 3);
  addFormFactor( 421,-325  ,2,-2, 4, 3);
  addFormFactor( 421,-10323,1,-2, 4, 3);
  addFormFactor( 421,-20323,1,-2, 4, 3);
  addFormFactor( 421,-10321,0,-2, 4, 3);
  // D+ to d dbar I=0
  addFormFactor( 411, 221  ,0,-1, 4, 1);
  addFormFactor( 411, 331  ,0,-1, 4, 1);
  addFormFactor( 411, 223  ,1,-1, 4, 1);
  addFormFactor( 411, 225  ,2,-1, 4, 1);
  addFormFactor( 411, 10223,1,-1, 4, 1);
  addFormFactor( 411, 20223,1,-1, 4, 1);
  addFormFactor( 411, 10221,0,-1, 4, 1);
  // D+ to d dbar I=1
  addFormFactor( 411, 111  ,0,-1, 4, 1);
  addFormFactor( 411, 113  ,1,-1, 4, 1);
  addFormFactor( 411, 115  ,2,-1, 4, 1);
  addFormFactor( 411, 10113,1,-1, 4, 1);
  addFormFactor( 411, 20113,1,-1, 4, 1);
  addFormFactor( 411, 10111,0,-1, 4, 1);
  // D+ to u dbar
  addFormFactor( 411, 211  ,0,-1, 4, 2);
  addFormFactor( 411, 213  ,1,-1, 4, 2);
  addFormFactor( 411, 215  ,2,-1, 4, 2);
  addFormFactor( 411, 10213,1,-1, 4, 2);
  addFormFactor( 411, 20213,1,-1, 4, 2);
  addFormFactor( 411, 10211,0,-1, 4, 2);
  // D+ to s dbar
  addFormFactor( 411,-311  ,0,-1, 4, 3);
  addFormFactor( 411,-313  ,1,-1, 4, 3);
  addFormFactor( 411,-315  ,2,-1, 4, 3);
  addFormFactor( 411,-10313,1,-1, 4, 3);
  addFormFactor( 411,-20313,1,-1, 4, 3);
  addFormFactor( 411,-10311,0,-1, 4, 3);
  // set the initial number of modes
  initialModes(numberOfFactors());
}			     

void ISGW2FormFactor::doinit() {
  ScalarFormFactor::doinit();
  // set up the quark masses
  _mquark.resize(5);
  _mquark[0]=_mdown;
  _mquark[1]=_mup;
  _mquark[2]=_mstrange;
  _mquark[3]=_mcharm;
  _mquark[4]=_mbottom;
  // value of alpha_S at the quark masses
  _alphaQ.resize(5);
  for(unsigned int ix=0;ix<5;++ix) {
    _alphaQ[ix]=alphaS(_mquark[ix],sqr(_mquark[ix]));
  }
  _beta1S0.resize(5,vector<Energy>(5));
  _mass1S0.resize(5,vector<Energy>(5));
  _beta3S1.resize(5,vector<Energy>(5));
  _beta1P .resize(5,vector<Energy>(5));
  _massPoh.resize(5,vector<Energy>(5));
  _massPth.resize(5,vector<Energy>(5));
  // set up the beta values
  _beta1S0[0][0] = _beta1S0ud;_beta3S1[0][0] = _beta3S1ud;_beta1P[0][0] = _beta1Pud;
  _beta1S0[1][0] = _beta1S0ud;_beta3S1[1][0] = _beta3S1ud;_beta1P[1][0] = _beta1Pud;
  _beta1S0[2][0] = _beta1S0us;_beta3S1[2][0] = _beta3S1us;_beta1P[2][0] = _beta1Pus;
  _beta1S0[3][0] = _beta1S0cu;_beta3S1[3][0] = _beta3S1cu;_beta1P[3][0] = _beta1Pcu;
  _beta1S0[4][0] = _beta1S0ub;_beta3S1[4][0] = _beta3S1ub;_beta1P[4][0] = _beta1Pub;
  _beta1S0[0][1] = _beta1S0ud;_beta3S1[0][1] = _beta3S1ud;_beta1P[0][1] = _beta1Pud;
  _beta1S0[1][1] = _beta1S0ud;_beta3S1[1][1] = _beta3S1ud;_beta1P[1][1] = _beta1Pud;
  _beta1S0[2][1] = _beta1S0us;_beta3S1[2][1] = _beta3S1us;_beta1P[2][1] = _beta1Pus;
  _beta1S0[3][1] = _beta1S0cu;_beta3S1[3][1] = _beta3S1cu;_beta1P[3][1] = _beta1Pcu;
  _beta1S0[4][1] = _beta1S0ub;_beta3S1[4][1] = _beta3S1ub;_beta1P[4][1] = _beta1Pub;
  _beta1S0[0][2] = _beta1S0us;_beta3S1[0][2] = _beta3S1us;_beta1P[0][2] = _beta1Pus;
  _beta1S0[1][2] = _beta1S0us;_beta3S1[1][2] = _beta3S1us;_beta1P[1][2] = _beta1Pus;
  _beta1S0[2][2] = _beta1S0ss;_beta3S1[2][2] = _beta3S1ss;_beta1P[2][2] = _beta1Pss;
  _beta1S0[3][2] = _beta1S0cs;_beta3S1[3][2] = _beta3S1cs;_beta1P[3][2] = _beta1Pcs;
  _beta1S0[4][2] = _beta1S0sb;_beta3S1[4][2] = _beta3S1sb;_beta1P[4][2] = _beta1Psb;
  _beta1S0[0][3] = _beta1S0cu;_beta3S1[0][3] = _beta3S1cu;_beta1P[0][3] = _beta1Pcu;
  _beta1S0[1][3] = _beta1S0cu;_beta3S1[1][3] = _beta3S1cu;_beta1P[1][3] = _beta1Pcu;
  _beta1S0[2][3] = _beta1S0cs;_beta3S1[2][3] = _beta3S1cs;_beta1P[2][3] = _beta1Pcs;
  _beta1S0[3][3] = _beta1S0cc;_beta3S1[3][3] = _beta3S1cc;_beta1P[3][3] = _beta1Pcc;
  _beta1S0[4][3] = _beta1S0bc;_beta3S1[4][3] = _beta3S1bc;_beta1P[4][3] = _beta1Pbc;
  _beta1S0[0][4] = _beta1S0ub;_beta3S1[0][4] = _beta3S1ub;_beta1P[0][4] = _beta1Pub;
  _beta1S0[1][4] = _beta1S0ub;_beta3S1[1][4] = _beta3S1ub;_beta1P[1][4] = _beta1Pub;
  _beta1S0[2][4] = _beta1S0sb;_beta3S1[2][4] = _beta3S1sb;_beta1P[2][4] = _beta1Psb;
  _beta1S0[3][4] = _beta1S0bc;_beta3S1[3][4] = _beta3S1bc;_beta1P[3][4] = _beta1Pbc;
  _beta1S0[4][4] = ZERO    ;_beta3S1[4][4] = ZERO    ;_beta1P[4][4] = ZERO   ;
  // set up the values of mbar
  // get the particle data objects
  tcPDPtr p1S0[5][5],p3S1[5][5];
  p1S0[0][0] = getParticleData(111); p3S1[0][0] = getParticleData(113);
  p1S0[1][0] = getParticleData(211); p3S1[1][0] = getParticleData(213);
  p1S0[2][0] = getParticleData(311); p3S1[2][0] = getParticleData(313);
  p1S0[3][0] = getParticleData(411); p3S1[3][0] = getParticleData(413);
  p1S0[4][0] = getParticleData(511); p3S1[4][0] = getParticleData(513);
  p1S0[0][1] = getParticleData(211); p3S1[0][1] = getParticleData(213);
  p1S0[1][1] = getParticleData(111); p3S1[1][1] = getParticleData(113);
  p1S0[2][1] = getParticleData(321); p3S1[2][1] = getParticleData(323);
  p1S0[3][1] = getParticleData(421); p3S1[3][1] = getParticleData(423);
  p1S0[4][1] = getParticleData(521); p3S1[4][1] = getParticleData(523);
  p1S0[0][2] = getParticleData(311); p3S1[0][2] = getParticleData(313);
  p1S0[1][2] = getParticleData(321); p3S1[1][2] = getParticleData(323);
  p1S0[2][2] = getParticleData(331); p3S1[2][2] = getParticleData(333);
  p1S0[3][2] = getParticleData(431); p3S1[3][2] = getParticleData(433);
  p1S0[4][2] = getParticleData(531); p3S1[4][2] = getParticleData(533);
  p1S0[0][3] = getParticleData(411); p3S1[0][3] = getParticleData(413);
  p1S0[1][3] = getParticleData(421); p3S1[1][3] = getParticleData(423);
  p1S0[2][3] = getParticleData(431); p3S1[2][3] = getParticleData(433);
  p1S0[3][3] = getParticleData(441); p3S1[3][3] = getParticleData(443);
  p1S0[4][3] = getParticleData(541); p3S1[4][3] = getParticleData(543);
  p1S0[0][4] = getParticleData(511); p3S1[0][4] = getParticleData(513);
  p1S0[1][4] = getParticleData(521); p3S1[1][4] = getParticleData(523);
  p1S0[2][4] = getParticleData(531); p3S1[2][4] = getParticleData(533);
  p1S0[3][4] = getParticleData(541); p3S1[3][4] = getParticleData(543);
  p1S0[4][4] = getParticleData(551); p3S1[4][4] = getParticleData(553);
  tcPDPtr p3P0[5][5],p3P1[5][5],p3P2[5][5],p1P1[5][5];
  p3P0[0][0] = getParticleData(10111); p3P1[0][0] = getParticleData(20113);
  p3P0[1][0] = getParticleData(10211); p3P1[1][0] = getParticleData(20213);
  p3P0[2][0] = getParticleData(10311); p3P1[2][0] = getParticleData(20313);
  p3P0[3][0] = getParticleData(10411); p3P1[3][0] = getParticleData(20413);
  p3P0[4][0] = getParticleData(10511); p3P1[4][0] = getParticleData(20513);
  p3P0[0][1] = getParticleData(10211); p3P1[0][1] = getParticleData(20213);
  p3P0[1][1] = getParticleData(10111); p3P1[1][1] = getParticleData(20113);
  p3P0[2][1] = getParticleData(10321); p3P1[2][1] = getParticleData(20323);
  p3P0[3][1] = getParticleData(10421); p3P1[3][1] = getParticleData(20423);
  p3P0[4][1] = getParticleData(10521); p3P1[4][1] = getParticleData(20523);
  p3P0[0][2] = getParticleData(10311); p3P1[0][2] = getParticleData(20313);
  p3P0[1][2] = getParticleData(10321); p3P1[1][2] = getParticleData(20323);
  p3P0[2][2] = getParticleData(10331); p3P1[2][2] = getParticleData(20333);
  p3P0[3][2] = getParticleData(10431); p3P1[3][2] = getParticleData(20433);
  p3P0[4][2] = getParticleData(10531); p3P1[4][2] = getParticleData(20533);
  p3P0[0][3] = getParticleData(10411); p3P1[0][3] = getParticleData(20413);
  p3P0[1][3] = getParticleData(10421); p3P1[1][3] = getParticleData(20423);
  p3P0[2][3] = getParticleData(10431); p3P1[2][3] = getParticleData(20433);
  p3P0[3][3] = getParticleData(10441); p3P1[3][3] = getParticleData(20443);
  p3P0[4][3] = getParticleData(10541); p3P1[4][3] = getParticleData(20543);
  p3P0[0][4] = getParticleData(10511); p3P1[0][4] = getParticleData(20513);
  p3P0[1][4] = getParticleData(10521); p3P1[1][4] = getParticleData(20523);
  p3P0[2][4] = getParticleData(10531); p3P1[2][4] = getParticleData(20533);
  p3P0[3][4] = getParticleData(10541); p3P1[3][4] = getParticleData(20543);
  p3P0[4][4] = getParticleData(10551); p3P1[4][4] = getParticleData(20553);
  p1P1[0][0]=getParticleData(10113); p3P2[0][0]=getParticleData(115);
  p1P1[1][0]=getParticleData(10213); p3P2[1][0]=getParticleData(215);
  p1P1[2][0]=getParticleData(10313); p3P2[2][0]=getParticleData(315);
  p1P1[3][0]=getParticleData(10413); p3P2[3][0]=getParticleData(415);
  p1P1[4][0]=getParticleData(10513); p3P2[4][0]=getParticleData(515);
  p1P1[0][1]=getParticleData(10213); p3P2[0][1]=getParticleData(215);
  p1P1[1][1]=getParticleData(10113); p3P2[1][1]=getParticleData(115);
  p1P1[2][1]=getParticleData(10323); p3P2[2][1]=getParticleData(325);
  p1P1[3][1]=getParticleData(10423); p3P2[3][1]=getParticleData(425);
  p1P1[4][1]=getParticleData(10523); p3P2[4][1]=getParticleData(525);
  p1P1[0][2]=getParticleData(10313); p3P2[0][2]=getParticleData(315);
  p1P1[1][2]=getParticleData(10323); p3P2[1][2]=getParticleData(325);
  p1P1[2][2]=getParticleData(10333); p3P2[2][2]=getParticleData(335);
  p1P1[3][2]=getParticleData(10433); p3P2[3][2]=getParticleData(435);
  p1P1[4][2]=getParticleData(10533); p3P2[4][2]=getParticleData(535);
  p1P1[0][3]=getParticleData(10413); p3P2[0][3]=getParticleData(415);
  p1P1[1][3]=getParticleData(10423); p3P2[1][3]=getParticleData(425);
  p1P1[2][3]=getParticleData(10433); p3P2[2][3]=getParticleData(435);
  p1P1[3][3]=getParticleData(10443); p3P2[3][3]=getParticleData(445);
  p1P1[4][3]=getParticleData(10543); p3P2[4][3]=getParticleData(545);
  p1P1[0][4]=getParticleData(10513); p3P2[0][4]=getParticleData(515);
  p1P1[1][4]=getParticleData(10523); p3P2[1][4]=getParticleData(525);
  p1P1[2][4]=getParticleData(10533); p3P2[2][4]=getParticleData(535);
  p1P1[3][4]=getParticleData(10543); p3P2[3][4]=getParticleData(545);
  p1P1[4][4]=getParticleData(10553); p3P2[4][4]=getParticleData(555);
  // calculate the masses
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      Energy m1S0,m3S1,m3P0,m3P1,m3P2,m1P1;
      if(!p1S0[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 1S0 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m1S0 = ZERO;
      }
      else {
	m1S0 = p1S0[ix][iy]->mass();
      }
      if(!p3S1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3S1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3S1 = ZERO;
      }
      else {
	m3S1 = p3S1[ix][iy]->mass();
      }
      if(!p3P0[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P0 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P0 = ZERO;
      }
      else {
	m3P0 = p3P0[ix][iy]->mass();
      }
      if(!p3P1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P1 = ZERO;
      }
      else {
	m3P1 = p3P1[ix][iy]->mass();
      }
      if(!p3P2[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P2 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P2 = ZERO;
      }
      else {
	m3P2 = p3P2[ix][iy]->mass();
      }
      if(!p1P1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 1P1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m1P1 = ZERO;
      }
      else {
	m1P1 = p1P1[ix][iy]->mass();
      }
      // 1S0
      _mass1S0[ix][iy] = 0.75 *m3S1+0.25 *m1S0;
      //  1p 1/2
      _massPoh[ix][iy] = 0.75 *m3P1+0.25 *m3P0;
      //  1p 3/2
      _massPth[ix][iy] = 0.625*m3P2+0.375*m1P1;
    }
  }
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Energy mtemp = (4.*_massPoh[ix][iy]+8.*_massPth[ix][iy])/12.;
      _massPoh[ix][iy]=mtemp;
      _massPth[ix][iy]=mtemp;
    }
  }
}

void ISGW2FormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mdown,GeV)  << ounit(_mup,GeV) << ounit(_mstrange,GeV) 
     << ounit(_mcharm,GeV) << ounit(_mbottom,GeV) << ounit(_beta1S0ud,GeV)
     << ounit(_beta1S0us,GeV) << ounit(_beta1S0ss,GeV) << ounit(_beta1S0cu,GeV) 
     << ounit(_beta1S0cs,GeV) << ounit(_beta1S0ub,GeV) << ounit(_beta1S0sb,GeV) 
     << ounit(_beta1S0cc,GeV) << ounit(_beta1S0bc,GeV) << ounit(_beta3S1ud,GeV) 
     << ounit(_beta3S1us,GeV) << ounit(_beta3S1ss,GeV) << ounit(_beta3S1cu,GeV) 
     << ounit(_beta3S1cs,GeV) << ounit(_beta3S1ub,GeV) << ounit(_beta3S1sb,GeV) 
     << ounit(_beta3S1cc,GeV) << ounit(_beta3S1bc,GeV) << ounit(_beta1Pud ,GeV) 
     << ounit(_beta1Pus ,GeV) << ounit(_beta1Pss ,GeV) << ounit(_beta1Pcu ,GeV) 
     << ounit(_beta1Pcs ,GeV) << ounit(_beta1Pub ,GeV) << ounit(_beta1Psb ,GeV) 
     << ounit(_beta1Pcc ,GeV) << ounit(_beta1Pbc ,GeV)
     << _alphamuQM  << _CfDrho << _CfDKstar << _CfDsKstar << _CfDsphi 
     << _CfBrho << _CfBDstar << _CfBsKstar << _CfBsDstar << _CfBcDstar << _CfBcpsi
     << _CfBcBsstar << _CfBcBstar << _thetaeta << ounit(_mquark,GeV) << _alphaQ 
     << ounit(_beta1S0,GeV) << ounit(_mass1S0,GeV) << ounit(_beta3S1,GeV) 
     << ounit(_beta1P,GeV) << ounit(_massPoh,GeV) << ounit(_massPth,GeV) << _includeaW;
}

void ISGW2FormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mdown,GeV) >> iunit(_mup,GeV) >> iunit(_mstrange,GeV) 
     >> iunit(_mcharm,GeV) >> iunit(_mbottom,GeV) >> iunit(_beta1S0ud,GeV) 
     >> iunit(_beta1S0us,GeV) >> iunit(_beta1S0ss,GeV) >> iunit(_beta1S0cu,GeV) 
     >> iunit(_beta1S0cs,GeV) >> iunit(_beta1S0ub,GeV) >> iunit(_beta1S0sb,GeV) 
     >> iunit(_beta1S0cc,GeV) >> iunit(_beta1S0bc,GeV) >> iunit(_beta3S1ud,GeV) 
     >> iunit(_beta3S1us,GeV) >> iunit(_beta3S1ss,GeV) >> iunit(_beta3S1cu,GeV) 
     >> iunit(_beta3S1cs,GeV) >> iunit(_beta3S1ub,GeV) >> iunit(_beta3S1sb,GeV) 
     >> iunit(_beta3S1cc,GeV) >> iunit(_beta3S1bc,GeV) >> iunit(_beta1Pud ,GeV) 
     >> iunit(_beta1Pus ,GeV) >> iunit(_beta1Pss ,GeV) >> iunit(_beta1Pcu ,GeV) 
     >> iunit(_beta1Pcs ,GeV) >> iunit(_beta1Pub ,GeV) >> iunit(_beta1Psb ,GeV) 
     >> iunit(_beta1Pcc ,GeV) >> iunit(_beta1Pbc ,GeV) 
     >> _alphamuQM >> _CfDrho >> _CfDKstar >> _CfDsKstar >> _CfDsphi 
     >> _CfBrho >> _CfBDstar >> _CfBsKstar >> _CfBsDstar >> _CfBcDstar >> _CfBcpsi
     >> _CfBcBsstar >> _CfBcBstar >> _thetaeta >> iunit(_mquark,GeV) >> _alphaQ 
     >> iunit(_beta1S0,GeV) >> iunit(_mass1S0,GeV) >> iunit(_beta3S1,GeV) 
     >> iunit(_beta1P,GeV) >> iunit(_massPoh,GeV) >> iunit(_massPth,GeV) >> _includeaW;
}

ClassDescription<ISGW2FormFactor> ISGW2FormFactor::initISGW2FormFactor;
// Definition of the static class description member.

void ISGW2FormFactor::Init() {

  static ClassDocumentation<ISGW2FormFactor> documentation
    ("The ISGW2FormFactor class implements the ISGW2 model of "
     "PRD52, 2783.",
     "The ISGW2 form factor model \\cite{Scora:1995ty} was used.",
     "\\bibitem{Scora:1995ty} D.~Scora and N.~Isgur,"
     "Phys.\\ Rev.\\  D {\\bf 52} (1995) 2783 [arXiv:hep-ph/9503486].\n"
     "%%CITATION = PHRVA,D52,2783;%%\n");

  static Parameter<ISGW2FormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mdown, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mup, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mstrange, GeV, 0.55*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mcharm, GeV, 1.82*GeV, ZERO, 3.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mbottom, GeV, 5.20*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ud
    ("Beta1S0ud",
     "The beta wavefunction parameter for the ud meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ud, GeV, 0.41*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0us
    ("Beta1S0us",
     "The beta wavefunction parameter for the us meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0us, GeV, 0.44*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ss
    ("Beta1S0ss",
     "The beta wavefunction parameter for the ss meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ss, GeV, 0.53*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cu
    ("Beta1S0cu",
     "The beta wavefunction parameter for the cu meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cu, GeV, 0.45*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cs
    ("Beta1S0cs",
     "The beta wavefunction parameter for the cs meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cs, GeV, 0.56*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ub
    ("Beta1S0ub",
     "The beta wavefunction parameter for the ub meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ub, GeV, 0.43*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0sb
    ("Beta1S0sb",
     "The beta wavefunction parameter for the sb meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0sb, GeV, 0.54*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cc
    ("Beta1S0cc",
     "The beta wavefunction parameter for the cc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cc, GeV, 0.88*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0bc
    ("Beta1S0bc",
     "The beta wavefunction parameter for the bc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0bc, GeV, 0.92*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pud
    ("Beta1Pud",
     "The beta wavefunction parameter for the ud meson in the 1P level",
     &ISGW2FormFactor::_beta1Pud, GeV, 0.28*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pus
    ("Beta1Pus",
     "The beta wavefunction parameter for the us meson in the 1P level",
     &ISGW2FormFactor::_beta1Pus, GeV, 0.30*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pss
    ("Beta1Pss",
     "The beta wavefunction parameter for the ss meson in the 1P level",
     &ISGW2FormFactor::_beta1Pss, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcu
    ("Beta1Pcu",
     "The beta wavefunction parameter for the cu meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcu, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcs
    ("Beta1Pcs",
     "The beta wavefunction parameter for the cs meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcs, GeV, 0.38*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pub
    ("Beta1Pub",
     "The beta wavefunction parameter for the ub meson in the 1P level",
     &ISGW2FormFactor::_beta1Pub, GeV, 0.35*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Psb
    ("Beta1Psb",
     "The beta wavefunction parameter for the sb meson in the 1P level",
     &ISGW2FormFactor::_beta1Psb, GeV, 0.41*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcc
    ("Beta1Pcc",
     "The beta wavefunction parameter for the cc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcc, GeV, 0.52*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pbc
    ("Beta1Pbc",
     "The beta wavefunction parameter for the bc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pbc, GeV, 0.60*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ud
    ("Beta3S1ud",
     "The beta wavefunction parameter for the ud meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ud, GeV, 0.30*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1us
    ("Beta3S1us",
     "The beta wavefunction parameter for the us meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1us, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ss
    ("Beta3S1ss",
     "The beta wavefunction parameter for the ss meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ss, GeV, 0.37*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cu
    ("Beta3S1cu",
     "The beta wavefunction parameter for the cu meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cu, GeV, 0.38*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cs
    ("Beta3S1cs",
     "The beta wavefunction parameter for the cs meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cs, GeV, 0.44*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ub
    ("Beta3S1ub",
     "The beta wavefunction parameter for the ub meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ub, GeV, 0.40*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1sb
    ("Beta3S1sb",
     "The beta wavefunction parameter for the sb meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1sb, GeV, 0.49*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cc
    ("Beta3S1cc",
     "The beta wavefunction parameter for the cc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cc, GeV, 0.62*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1bc
    ("Beta3S1bc",
     "The beta wavefunction parameter for the bc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1bc, GeV, 0.75*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceAlphaCutOff
    ("AlphaCutOff",
     "The value of the strong coupling constnat at the cut-off",
     &ISGW2FormFactor::_alphamuQM, 0.6, 0.0, 10.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDrho
    ("CfDrho",
     "The relativistic correction factor for D -> rho",
     &ISGW2FormFactor::_CfDrho, 0.889, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDKstar
    ("CfDKstar",
     "The relativistic correction factor for D -> Kstar",
     &ISGW2FormFactor::_CfDKstar, 0.928, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsKstar
    ("CfDsKstar",
     "The relativistic correction factor for Ds -> Kstar",
     &ISGW2FormFactor::_CfDsKstar, 0.873, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsphi
    ("CfDsphi",
     "The relativistic correction factor for Ds -> phi",
     &ISGW2FormFactor::_CfDsphi, 0.911, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBrho
    ("CfBrho",
     "The relativistic correction factor for B -> rho",
     &ISGW2FormFactor::_CfBrho, 0.905, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBDstar
    ("CfBDstar",
     "The relativistic correction factor for B -> Dstar",
     &ISGW2FormFactor::_CfBDstar, 0.989, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsKstar
    ("CfBsKstar",
     "The relativistic correction factor for Bs -> Kstar",
     &ISGW2FormFactor::_CfBsKstar, 0.892, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsDstar
    ("CfBsDstar",
     "The relativistic correction factor for Bs -> Dstar",
     &ISGW2FormFactor::_CfBsDstar, 0.984, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcDstar
    ("CfBcDstar",
     "The relativistic correction factor for Bc -> Dstar",
     &ISGW2FormFactor::_CfBcDstar, 0.868, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcpsi
    ("CfBcpsi",
     "The relativistic correction factor for Bc -> psi",
     &ISGW2FormFactor::_CfBcpsi, 0.967, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBsstar
    ("CfBcBsstar",
     "The relativistic correction factor for Bc -> Bsstar",
     &ISGW2FormFactor::_CfBcBsstar, 1.000, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBstar
    ("CfBcBstar",
     "The relativistic correction factor for Bc -> Bstar",
     &ISGW2FormFactor::_CfBcBstar, 1.000, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ISGW2FormFactor::_thetaeta, -Constants::pi/9., -Constants::pi, Constants::pi,
     false, false, true);

  static Switch<ISGW2FormFactor,bool> interfaceIncludeaW
    ("IncludeaW",
     "Include the a(omega) piece of the Cji factor",
     &ISGW2FormFactor::_includeaW, true, false, false);
  static SwitchOption interfaceIncludeaWInclude
    (interfaceIncludeaW,
     "Yes",
     "Include the factor",
     true);
  static SwitchOption interfaceIncludeaWDoNot
    (interfaceIncludeaW,
     "No",
     "Do not include the factor",
     false);

}

// member which does the work
void ISGW2FormFactor::formFactor(Energy2 q2, unsigned int iloc, int, int id1,
				 Energy mY,
				 Energy mX, Complex & f1,Complex & f2,Complex & f3,
				 Complex & f4) const {
  useMe();
  // get the flavours of the quarks etc
  int jspin,spect,inquark,outquark;
  formFactorInfo(iloc,jspin,spect,inquark,outquark);
  int ifl0(abs(inquark)),ifl1(abs(outquark)),ifls(abs(spect));
  // determine the multiplet
  int ispin(abs(id1)/1000);
  // masses of the quarks
  Energy mQ(_mquark[ifl0-1]),mq(_mquark[ifl1-1]),ms(_mquark[ifls-1]);
  // of the mesons
  Energy mtildeX(mq+ms),mtildeY(mQ+ms),mup(mq*mQ/(mQ+mq)),mum(mq*mQ/(mQ-mq));
  // wavefunction parameters for the mesons
  Energy betaX(ZERO),mbarX(ZERO),
    betaY(_beta1S0[ifl0-1][ifls-1]),mbarY(_mass1S0[ifl0-1][ifls-1]);
  double Cf(1.);
  // the wavefunction parameter for the outgoing meson
  // 1S0
  if(ispin==0&&jspin==0) {
    betaX=_beta1S0[ifl1-1][ifls-1];
    mbarX=_mass1S0[ifl1-1][ifls-1];
  }
  // 3S1
  else if(ispin==0&&jspin==1) {
    betaX = _beta3S1[ifl1-1][ifls-1];
    mbarX = _mass1S0[ifl1-1][ifls-1];
    // set the relativistic correction parameter
    // decaying b
    if(ifl0==5) {
      if(ifls<3)       Cf = ifl1<3  ? _CfBrho    : _CfBDstar;
      else if(ifls==3) Cf = ifl1==4 ? _CfBsDstar : _CfBsKstar;
      else if(ifls==4) Cf = ifl1==4 ? _CfBcpsi   : _CfBcDstar;
    }
    // decaying D
    else if(ifl0==4) {
      if(ifls<3)       Cf = ifl1<3  ? _CfDrho    : _CfDKstar;
      else if(ifls==3) Cf = ifl1<3  ? _CfDsKstar : _CfDsphi;
      else if(ifls==5) Cf = ifl1<3  ? _CfBcBstar : _CfBcBsstar;
    } 
  }
  else if(ispin==10&&jspin==0) {
    betaX=_beta1P[ifl1-1][ifls-1];
    mbarX=_massPoh[ifl1-1][ifls-1];
  }
  // 1 3/2 P 1 (1 P1)
  else if((ispin==0&&jspin==2)||(ispin==10&&jspin==1)) {
    betaX = _beta1P[ifl1-1][ifls-1];
    mbarX=_massPth[ifl1-1][ifls-1];
  }
  // 1 1/2 P1 ( 3 P1) 
  else if(ispin==20&&jspin==1) {
    betaX = _beta1P[ifl1-1][ifls-1];
    mbarX=_massPoh[ifl1-1][ifls-1];
  }
  else {
    throw Exception() << "ISGWS2FormFactor::formFactor" 
		      << " unknown multiplet" << Exception::abortnow;
  }
  Energy2 beta2XY(0.5*(betaX*betaX+betaY*betaY));
  // number of active flavours
  int Nf  = ifl0-1;
  int Nfp = ifl1==2 ? 0 : ifl1-1;
  // first piece of the f_n function
  double betar(betaX*betaY/beta2XY),fn(sqrt(mtildeX/mtildeY)*betar*sqrt(betar));
  // q dependent piece
  Energy2 tm((mY-mX)*(mY-mX)),tmmt(tm-q2);
  // radius parameter
  InvEnergy2 r2(0.75/mQ/mq+1.5*ms*ms/mbarX/mbarY/beta2XY
		+16./mbarX/mbarY/(33.-2.*Nfp)*log(_alphamuQM/_alphaQ[ifl1-1]));
  // the parameters for the form-factor depenedent piece
  double rmbmtY(sqrt(mbarY/mtildeY)),rmbmtX(sqrt(mbarX/mtildeX));
  // work out wtilde
  double wt(1.+0.5*tmmt/mbarX/mbarY);
  // storage of the form factors
  Energy f(ZERO);
  InvEnergy g(ZERO),appam(ZERO),apmam(ZERO);
  InvEnergy2 h(ZERO),bp(ZERO),bm(ZERO);
  double fpmfm(0.),fppfm(0.),k(0.);
  // scalar and vector from 1S levels
  if(ispin==0) {
    // parameters for the beta functions
    double asopi(alphaS(mq,mq*mQ)/Constants::pi),w(1.+0.5*tmmt/mX/mY);
    double aI(-6./(33.-2.*Nf)),rw(1./sqrt(w*w-1)*log(w+sqrt(w*w-1.)));
    double aLw(8./(33.-2.*Nfp)*(w*rw-1.)); 
    double cji(pow(_alphaQ[ifl0-1]/_alphaQ[ifl1-1],aI));
    if(_includeaW) cji*=pow(_alphaQ[ifl1-1]/_alphamuQM,aLw);
    double zji(mq/mQ); 
    double gamji(-2.*zji/(1.-zji)*log(zji)-2.),chiji(-1.-gamji/(1.-zji));
    // scalar
    if(jspin==0) {
      double fact((1.+1./12.*r2*tmmt));
      fn/=(fact*fact);
      fact = (1.-0.5*ms*mq/mup/mtildeX*betaY*betaY/beta2XY);
      fppfm = fn*rmbmtX/rmbmtY*cji*(1.+asopi*(gamji-2./3.*chiji))*
	(2.-mtildeX/mq*fact);
      fpmfm = fn*rmbmtY/rmbmtX*cji*(1.+asopi*(gamji+2./3.*chiji))*mtildeY/mq*fact;
    }
    else if(jspin==1) {
      // factors for the F and R functions
      double fact((1.+1./12.*r2*tmmt));
      fn/=(fact*fact);
      double betaapmam=1./3.-4./3./(1-zji)-chiji
	+gamji*(1.-2./3.*(1.+zji)/(1.-zji)/(1.-zji));
      double ftemp  = Cf*fn*rmbmtX*rmbmtY*cji*(1.+asopi*(-2./3.+gamji));
      double gtemp  = fn/rmbmtX/rmbmtY*cji*(1.+asopi*( 2./3.+gamji));
      double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY)*cji;
      double amtemp = fn/rmbmtX/rmbmtY*cji*(1.+asopi*betaapmam);
      // rest of the calculation
      f     =    ftemp*mtildeY*(1.+wt+0.5*ms*(wt-1.)/mup);
      g     =0.5*gtemp*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
      appam = aptemp*(ms*betaX*betaX/(1.+wt)/mq/mQ/beta2XY*
		      (1.-0.5*ms*betaX*betaX/mtildeY/beta2XY)
		      +asopi/mtildeY*(-1.-chiji+4./3./(1.-zji)
				      +2./3.*(1.+zji)/sqr(1.-zji)*gamji));
      apmam =-amtemp/mtildeX*(mtildeY/mQ
			      -0.5*ms*betaX*betaX/mup/beta2XY
			      +wt*ms*mtildeY*betaX*betaX/(wt+1.)/mq/mQ/beta2XY*
			      (1.-0.5*ms/mtildeY*betaX*betaX/beta2XY)); 
    }
    else if(jspin==2) {
      // factors for the F function
      double fact((1.+1./18.*r2*tmmt));
      fn*=betar/(fact*fact*fact);
      double htemp = fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
      double ktemp = fn*rmbmtX/rmbmtY;
      double bptemp(fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY*rmbmtY*rmbmtY));
      double bmtemp(fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY));
      // functions themselves
      double or2(sqrt(0.5));
      h = 0.5*htemp*ms*or2/mtildeY/betaY*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
      k = or2*ktemp*ms/betaY*(1.+wt);
      InvEnergy2 bppbm = 0.25*bptemp*or2*ms*ms/mq/mQ/mtildeY/betaY*betaX*betaX/beta2XY*
	(1.-0.5*ms/mtildeY*betaX*betaX/beta2XY);
      InvEnergy2 bpmbm = -or2*bmtemp*ms/mQ/mtildeX/betaY*
	(1.-0.5*ms*mQ/mup/mtildeY*betaX*betaX/beta2XY
	 +0.25*ms/mq*betaX*betaX/beta2XY*(1.-0.5*ms/mtildeY*betaX*betaX/beta2XY));
      // conversion
      bp = 0.5*(bppbm+bpmbm);
      bm = 0.5*(bppbm-bpmbm);
    }
  }
  // 1 3P0
  else if(ispin==10&&jspin==0) {
    fn*=betar;
    double fact=(1.+1./18.*r2*tmmt);
    fn/=(fact*fact*fact);
    fn *= sqrt(2./3.)*ms/betaY;
    fppfm =-fn*rmbmtX/rmbmtY;
    fpmfm = fn*rmbmtY/rmbmtX*mtildeY/mtildeX;
  }
  // 1 3/2 P1 ( 1 P1) 
  else if(ispin==10&&jspin==1) {
    // factors for the F and R functions
    double fact=(1.+1./18.*r2*tmmt);
    fn*=betar/(fact*fact*fact);
    double ftemp  = fn*rmbmtX*rmbmtY;
    double gtemp  = fn/rmbmtX/rmbmtY;
    double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
    double amtemp = fn/rmbmtX/rmbmtY;
    // light meson or onium
    if((ifls<3&&ifl1<3)||(ifls==ifl1)) {
      double oor2(sqrt(0.5));
      f     = oor2*ftemp*mtildeY*betaY*(1./mup
				       +ms*mtildeX/3./mq/betaY/betaY*(wt-1.)*(wt-1.));
      g     = oor2*gtemp *(0.25*mtildeY*betaY/mQ/mq/mtildeX+(wt-1.)*ms/6./mtildeX/betaY);
      appam = oor2*aptemp*ms/mtildeY/betaY*(1.-ms/mq+0.5*ms/mup*betaY*betaY/beta2XY);
      apmam = oor2*amtemp*ms/mq/betaY*((4.-wt)/3.
				       -0.5*ms*mq/mtildeX/mup*betaY*betaY/beta2XY);
    }
    // heavy meson
    else {
      double oor3(1./sqrt(3.));
      f     =-2.*ftemp*oor3*mtildeY*betaY*
	(1./mq+0.5*mtildeX*ms*(wt-1.)/betaY/betaY*
	 (0.5*(wt+1.)/mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY)); 
      g     =-0.5*gtemp*oor3*(0.5*(1.+wt)+0.5*betaY*betaY*mtildeY/ms/mq/mQ)*ms
	/betaY/mtildeX;
      appam =-0.5*aptemp/oor3*ms/betaY/mtildeY*
	(1.-ms/3./mq-ms/3.*betaY*betaY/beta2XY*
	 (0.5/mum-1./mup));
      apmam =-0.5*amtemp*oor3*ms/betaY/mtildeX*
	((2.-wt)*mtildeX/mq+ms*betaY*betaY/beta2XY*(0.5/mum-1./mup));
    }
  }
  // 1 1/2 P 1 (3 P1)
  else if(ispin==20&&jspin==1) {
    // factors for the F and R functions
    double fact=(1.+1./18.*r2*tmmt);
    fn*=betar/(fact*fact*fact);
    double ftemp  = fn*rmbmtX*rmbmtY;
    double gtemp  = fn/rmbmtX/rmbmtY;
    double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
    double amtemp = fn/rmbmtX/rmbmtY;
    // light meson
    if( ( ifls<3 && ifl1<3 ) || ( ifl1==ifls ) ) {
      f     =-ftemp *mtildeY*betaY*(1./mum
				   +ms*mtildeX*(wt-1.)/betaY/betaY*
				   ((5.+wt)/6./mq-0.5/mum*ms/mtildeX*
				    betaY*betaY/beta2XY)); 
      g     =-gtemp *0.5*ms/mtildeX/betaY*(5.+wt)/6.;
      appam =-aptemp*0.5*ms*mtildeX/mq/mtildeY/betaY*
	(1.-0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY); 
      apmam = amtemp*0.5*ms/mq/betaY*((wt+2.)/3.
				      -0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY);
    }
    // heavy meson
    else {
      double r2o3(sqrt(2./3.));
      f     =  ftemp*r2o3*mtildeY*betaY*(0.5/mq-1.5/mQ+ms*mtildeX*(wt-1.)/betaY/betaY*
					 (1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY));
      g     =  gtemp *0.5*r2o3*ms/betaY/mtildeX*(1.-0.25*betaY*betaY*mtildeY/ms/mq/mQ);
      appam =  aptemp*0.5*r2o3*ms*ms*betaX*betaX/mtildeY/mq/betaY/beta2XY; 
      apmam = -amtemp*r2o3*ms/mtildeX/betaY*(1.+0.5*ms*betaX*betaX/mq/beta2XY);
    }
  }
  else {
    throw Exception() << "ISGWS2FormFactor::formFactor" 
		      << " unknown multiplet" << Exception::abortnow;
  }
  // the final manipulations
  if(jspin==0) {
    double fp,fm;
    fp = 0.5*(fppfm+fpmfm);
    fm = 0.5*(fppfm-fpmfm);
    // convert to the standard form
    f1 = q2/(mY*mY-mX*mX)*fm+fp;
    f2 = fp;
  }
  else if(jspin==1) {
    InvEnergy ap(0.5*(appam+apmam)),am(0.5*(appam-apmam));
    // convert to the standard notation
    Energy msum(mX+mY),mdiff(mY-mX);
    Complex ii(0.,1.);
    f2 = -ii*f/msum;
    f3 = +Complex(ii*ap*msum);
    f1 = -Complex(ii*0.5/mX*(am*q2+ii*msum*f2-ii*mdiff*f3));
    f4 =  Complex(ii*g*msum);
  }
  else if(jspin==2) {
    Energy msum(mX+mY);
    f1 = h*sqr(msum);
    f2 = k;
    f3 = bp*sqr(msum);
    f4 = bm*sqr(msum);
  }
  // special for mixing
  double fact;
  if(id1==ParticleID::eta) {
    if(ifl1==3&&ifls==3){fact=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
    else{fact=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
  else if(id1==ParticleID::etaprime) {
    if(ifl1==3&&ifls==3){fact=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
    else{fact=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
  else if(ifl1==ifls&&ifl1<3) {
    if(abs(ifl1)==1&&int(id1/10)%100==1){fact=-sqrt(0.5);}
    else{fact=sqrt(0.5);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
}

// form-factor for scalar to scalar
void ISGW2FormFactor::ScalarScalarFormFactor(Energy2 q2, unsigned int iloc,int id0,
					     int id1, Energy mY, Energy mX,
					     Complex & f0,Complex & fp) const {
  Complex d1(0.),d2(0.);
  formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
}

// form-factor for scalar to vector
void ISGW2FormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, 
					     int id1, Energy mY, Energy mX,
					     Complex & A0,Complex & A1,
					     Complex & A2,Complex & V) const {
  formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);
}


// form-factor for scalar to tensor
void ISGW2FormFactor::
ScalarTensorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1, 
		       Energy mY, Energy mX, complex<InvEnergy2> & h,
		       Complex & k,complex<InvEnergy2> & bp,
		       complex<InvEnergy2> & bm) const {
  Complex f1,f2,f3,f4;
  formFactor(q2,iloc,id0,id1,mY,mX,f1,f2,f3,f4);
  Energy msum(mX+mY);
  h = f1/sqr(msum);
  k = f2;
  bp = f3/sqr(msum);
  bm = f4/sqr(msum);
}

void ISGW2FormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ISGW2FormFactor " << name() << "\n";
  output << "newdef " << name() << ":DownMass "    << _mdown/GeV    << "\n";
  output << "newdef " << name() << ":UpMass "      << _mup/GeV      << "\n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << "\n";
  output << "newdef " << name() << ":CharmMass "   << _mcharm/GeV   << "\n";
  output << "newdef " << name() << ":BottomMass "  << _mbottom/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ud " << _beta1S0ud/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0us " << _beta1S0us/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ss " << _beta1S0ss/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cu " << _beta1S0cu/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cs " << _beta1S0cs/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ub " << _beta1S0ub/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0sb " << _beta1S0sb/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cc " << _beta1S0cc/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0bc " << _beta1S0bc/GeV  << "\n";
  output << "newdef " << name() << ":Beta1Pud  " << _beta1Pud/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pus  " << _beta1Pus/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pss  " << _beta1Pss/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcu  " << _beta1Pcu/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcs  " << _beta1Pcs/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pub  " << _beta1Pub/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Psb  " << _beta1Psb/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcc  " << _beta1Pcc/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pbc  " << _beta1Pbc/GeV   << "\n";
  output << "newdef " << name() << ":Beta3S1ud " << _beta3S1ud/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1us " << _beta3S1us/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1ss " << _beta3S1ss/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cu " << _beta3S1cu/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cs " << _beta3S1cs/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1ub " << _beta3S1ub/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1sb " << _beta3S1sb/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cc " << _beta3S1cc/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1bc " << _beta3S1bc/GeV  << "\n";
  output << "newdef " << name() << ":AlphaCutOff " << _alphamuQM    << "\n";
  output << "newdef " << name() << ":CfDrho "      << _CfDrho       << "\n";
  output << "newdef " << name() << ":CfDKstar "    << _CfDKstar     << "\n";
  output << "newdef " << name() << ":CfDsKstar "   << _CfDsKstar    << "\n";
  output << "newdef " << name() << ":CfDsphi "     << _CfDsphi      << "\n";
  output << "newdef " << name() << ":CfBrho "      << _CfBrho       << "\n";
  output << "newdef " << name() << ":CfBDstar "    << _CfBDstar     << "\n";
  output << "newdef " << name() << ":CfBsKstar "   << _CfBsKstar    << "\n";
  output << "newdef " << name() << ":CfBsDstar "   << _CfBsDstar    << "\n";
  output << "newdef " << name() << ":CfBcDstar "   << _CfBcDstar    << "\n";
  output << "newdef " << name() << ":CfBcpsi "     << _CfBcpsi      << "\n";
  output << "newdef " << name() << ":CfBcBsstar "  << _CfBcBsstar   << "\n";
  output << "newdef " << name() << ":CfBcBstar "   << _CfBcBstar    << "\n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
