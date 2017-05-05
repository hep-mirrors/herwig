// -*- C++ -*-
//
// MelikhovStechFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovStechFormFactor class.
//

#include "MelikhovStechFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

MelikhovStechFormFactor::MelikhovStechFormFactor() 
: _fplus0(42,0.),_sigma1fp(42,0.),_sigma2fp(42,0.),
  _f00   (42,0.),_sigma1f0(42,0.),_sigma2f0(42,0.),
  _fT0   (42,0.),_sigma1fT(42,0.),_sigma2fT(42,0.),
  _V0    (42,0.),_sigma1V0(42,0.),_sigma2V0(42,0.),
  _A00   (42,0.),_sigma1A0(42,0.),_sigma2A0(42,0.),
  _A10   (42,0.),_sigma1A1(42,0.),_sigma2A1(42,0.),
  _A20   (42,0.),_sigma1A2(42,0.),_sigma2A2(42,0.),
  _T10   (42,0.),_sigma1T1(42,0.),_sigma2T1(42,0.),
  _T20   (42,0.),_sigma1T2(42,0.),_sigma2T2(42,0.),
  _T30   (42,0.),_sigma1T3(42,0.),_sigma2T3(42,0.),
  _massP(42,ZERO), _massV(42,ZERO) {
  // form factors for D to K
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  for(unsigned int ix=0;ix<2;++ix) {
    _fplus0[ix] = 0.78    ; _sigma1fp[ix] = 0.24    ;_sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.78    ; _sigma1f0[ix] = 0.38    ;_sigma2f0[ix] = 0.46;
    _fT0[ix]    = 0.75    ; _sigma1fT[ix] = 0.27    ;_sigma2fT[ix] = 0.00;
    _massP[ix]  = 1.97*GeV; _massV[ix]    = 2.11*GeV;
  }
  // form factors for D to K*
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  for(unsigned int ix=2;ix<4;++ix) {
    _V0[ix]    = 1.03    ; _sigma1V0[ix] = 0.27    ; _sigma2V0[ix] = 0.00;
    _A00[ix]   = 0.76    ; _sigma1A0[ix] = 0.17    ; _sigma2A0[ix] = 0.00;
    _A10[ix]   = 0.66    ; _sigma1A1[ix] = 0.30    ; _sigma2A1[ix] = 0.20;
    _A20[ix]   = 0.49    ; _sigma1A2[ix] = 0.67    ; _sigma2A2[ix] = 0.16;
    _T10[ix]   = 0.78    ; _sigma1T1[ix] = 0.25    ; _sigma2T1[ix] = 0.00;
    _T20[ix]   = 0.78    ; _sigma1T2[ix] = 0.02    ; _sigma2T2[ix] = 1.80;
    _T30[ix]   = 0.45    ; _sigma1T3[ix] = 1.23    ; _sigma2T3[ix] = 0.34;
    _massP[ix] = 1.97*GeV; _massV[ix]    = 2.11*GeV;
  }
  // form factors for D to pi
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  for(unsigned int ix=4;ix<8;++ix) {
    _fplus0[ix] = 0.69    ; _sigma1fp[ix] = 0.30; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.69    ; _sigma1f0[ix] = 0.54; _sigma2f0[ix] = 0.32;
    _fT0[ix]    = 0.60    ; _sigma1fT[ix] = 0.34; _sigma2fT[ix] = 0.00;
    _massP[ix]  = 1.87*GeV; _massV[ix]    = 2.01*GeV;
  }
  // form factors for D to rho
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  for(unsigned int ix=8;ix<12;++ix) {
    _V0[ix]  = 0.90; _sigma1V0[ix] = 0.46; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.66; _sigma1A0[ix] = 0.36; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.59; _sigma1A1[ix] = 0.50; _sigma2A1[ix] = 0.00;
    _A20[ix] = 0.49; _sigma1A2[ix] = 0.89; _sigma2A2[ix] = 0.00;
    _T10[ix] = 0.66; _sigma1T1[ix] = 0.44; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.66; _sigma1T2[ix] = 0.38; _sigma2T2[ix] = 0.50;
    _T30[ix] = 0.31; _sigma1T3[ix] = 1.10; _sigma2T3[ix] = 0.17;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // B to D
  addFormFactor(-521,421,0,-2,5,4);
  addFormFactor(-511,411,0,-2,5,4);
  for(unsigned int ix=12;ix<14;++ix) {
    _fplus0[ix] = 0.67; _sigma1fp[ix] = 0.57; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.67; _sigma1f0[ix] = 0.78; _sigma2f0[ix] = 0.00;
    _fT0[ix]    = 0.69; _sigma1fT[ix] = 0.56; _sigma2fT[ix] = 0.00;
    _massP[ix]  = 6.4*GeV; _massV[ix] = 6.4*GeV;
  }
  // B to D*
  addFormFactor(-521,423,1,-2,5,4);
  addFormFactor(-511,413,1,-1,5,4);
  for(unsigned int ix=14;ix<16;++ix) {
    _V0[ix]  = 0.76; _sigma1V0[ix] = 0.57; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.69; _sigma1A0[ix] = 0.58; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.66; _sigma1A1[ix] = 0.78; _sigma2A1[ix] = 0.00;
    _A20[ix] = 0.62; _sigma1A2[ix] = 1.40; _sigma2A2[ix] = 0.41;
    _T10[ix] = 0.68; _sigma1T1[ix] = 0.57; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.68; _sigma1T2[ix] = 0.64; _sigma2T2[ix] = 0.00;
    _T30[ix] = 0.33; _sigma1T3[ix] = 1.46; _sigma2T3[ix] = 0.50;
    _massP[ix] = 6.4*GeV; _massV[ix] = 6.4*GeV;
  }
  // B to K
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-1,5,3);
  for(unsigned int ix=16;ix<18;++ix) {
    _fplus0[ix] = 0.36; _sigma1fp[ix] = 0.43; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.36; _sigma1f0[ix] = 0.70; _sigma2f0[ix] = 0.27;
    _fT0[ix]    = 0.35; _sigma1fT[ix] = 0.43; _sigma2fT[ix] = 0.00;
    _massP[ix] = 5.37*GeV; _massV[ix] = 5.42*GeV;
  }
  // B to K*
  addFormFactor(-521,-323,1,-2,5,3);
  addFormFactor(-511,-313,1,-1,5,3);
  for(unsigned int ix=18;ix<20;++ix) {
    _V0[ix]  = 0.44; _sigma1V0[ix] = 0.45; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.45; _sigma1A0[ix] = 0.46; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.36; _sigma1A1[ix] = 0.64; _sigma2A1[ix] = 0.36;
    _A20[ix] = 0.32; _sigma1A2[ix] = 1.23; _sigma2A2[ix] = 0.38;
    _T10[ix] = 0.39; _sigma1T1[ix] = 0.45; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.39; _sigma1T2[ix] = 0.72; _sigma2T2[ix] = 0.62;
    _T30[ix] = 0.27; _sigma1T3[ix] = 1.31; _sigma2T3[ix] = 0.41;
    _massP[ix] = 5.37*GeV; _massV[ix] = 5.42*GeV;
  }
  // B to pi
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 211,0,-1,5,2);
  addFormFactor(-521,-211,0,-2,5,1);
  addFormFactor(-511, 111,0,-1,5,1);
  for(unsigned int ix=20;ix<24;++ix) {
    _fplus0[ix] = 0.29; _sigma1fp[ix] = 0.48; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.29; _sigma1f0[ix] = 0.76; _sigma2f0[ix] = 0.28;
    _fT0[ix]    = 0.28; _sigma1fT[ix] = 0.48; _sigma2fT[ix] = 0.00;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 213,1,-1,5,2);
  addFormFactor(-521,-213,1,-2,5,1);
  addFormFactor(-511, 113,1,-1,5,1);
  for(unsigned int ix=24;ix<28;++ix) {
    _V0[ix]  = 0.31; _sigma1V0[ix] = 0.59; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.30; _sigma1A0[ix] = 0.54; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.26; _sigma1A1[ix] = 0.73; _sigma2A1[ix] = 0.10;
    _A20[ix] = 0.24; _sigma1A2[ix] = 1.40; _sigma2A2[ix] = 0.50;
    _T10[ix] = 0.27; _sigma1T1[ix] = 0.60; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.27; _sigma1T2[ix] = 0.74; _sigma2T2[ix] = 0.19;
    _T30[ix] = 0.19; _sigma1T3[ix] = 1.42; _sigma2T3[ix] = 0.51;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // D_s to K
  addFormFactor(431,311,0,-3,4,1);
  addFormFactor(431,321,0,-3,4,2);
  for(unsigned int ix=28;ix<30;++ix) {
    _fplus0[ix] = 0.72; _sigma1fp[ix] = 0.20; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.72; _sigma1f0[ix] = 0.41; _sigma2f0[ix] = 0.70;
    _fT0[ix]    = 0.77; _sigma1fT[ix] = 0.24; _sigma2fT[ix] = 0.00;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // D_s to K*
  addFormFactor(431,313,1,-3,4,1);
  addFormFactor(431,323,1,-3,4,2);
  for(unsigned int ix=30;ix<32;++ix) {
    _V0[ix]  = 1.04; _sigma1V0[ix] =  0.24; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.67; _sigma1A0[ix] =  0.20; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.57; _sigma1A1[ix] =  0.29; _sigma2A1[ix] = 0.42;
    _A20[ix] = 0.42; _sigma1A2[ix] =  0.58; _sigma2A2[ix] = 0.00;
    _T10[ix] = 0.71; _sigma1T1[ix] =  0.22; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.71; _sigma1T2[ix] = -0.06; _sigma2T2[ix] = 0.44;
    _T30[ix] = 0.45; _sigma1T3[ix] =  1.08; _sigma2T3[ix] = 0.68;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // D_s to eta
  addFormFactor(431,221,0,-3,4,3);
  _fplus0[32] = 0.78; _sigma1fp[32] = 0.23; _sigma2fp[32] = 0.00;
  _f00[32]    = 0.78; _sigma1f0[32] = 0.33; _sigma2f0[32] = 0.38;
  _fT0[32]    = 0.80; _sigma1fT[32] = 0.24; _sigma2fT[32] = 0.00;
  _massP[32] = 1.97*GeV; _massV[32] = 2.11*GeV;
  // D_s to eta'
  addFormFactor(431,331,0,-3,4,3);
  _fplus0[33] = 0.78; _sigma1fp[33] = 0.23; _sigma2fp[33] = 0.00;
  _f00[33]    = 0.78; _sigma1f0[33] = 0.21; _sigma2f0[33] = 0.76;
  _fT0[33]    = 0.94; _sigma1fT[33] = 0.24; _sigma2fT[33] = 0.00;
  _massP[33] = 1.97*GeV; _massV[33] = 2.11*GeV;
  // D_s to phi
  addFormFactor(431,333,1,-3,4,3);
  _V0[34]  = 1.10; _sigma1V0[34] = 0.26; _sigma2V0[34] = 0.00;
  _A00[34] = 0.73; _sigma1A0[34] = 0.10; _sigma2A0[34] = 0.00;
  _A10[34] = 0.64; _sigma1A1[34] = 0.29; _sigma2A1[34] = 0.00;
  _A20[34] = 0.47; _sigma1A2[34] = 0.63; _sigma2A2[34] = 0.00;
  _T10[34] = 0.77; _sigma1T1[34] = 0.25; _sigma2T1[34] = 0.00;
  _T20[34] = 0.77; _sigma1T2[34] = 0.02; _sigma2T2[34] = 2.01;
  _T30[34] = 0.46; _sigma1T3[34] = 1.34; _sigma2T3[34] = 0.45;
  _massP[34] = 1.97*GeV; _massV[34] = 2.11*GeV;
  // B_s to K
  addFormFactor(-531, 311  ,0,-3,5,1);
  addFormFactor(-531, 321  ,0,-3,5,2);
  for(unsigned int ix=35;ix<37;++ix) {
    _fplus0[ix] = 0.31; _sigma1fp[ix] = 0.63; _sigma2fp[ix] = 0.33;
    _f00[ix]    = 0.31; _sigma1f0[ix] = 0.93; _sigma2f0[ix] = 0.70;
    _fT0[ix]    = 0.31; _sigma1fT[ix] = 0.61; _sigma2fT[ix] = 0.30;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B_s to K*
  addFormFactor(-531, 313  ,1,-3,5,1);
  addFormFactor(-531, 323  ,1,-3,5,2);
  for(unsigned int ix=37;ix<39;++ix) {
    _V0[ix]  = 0.38; _sigma1V0[ix] = 0.66; _sigma2V0[ix] = 0.30;
    _A00[ix] = 0.37; _sigma1A0[ix] = 0.60; _sigma2A0[ix] = 0.16;
    _A10[ix] = 0.29; _sigma1A1[ix] = 0.86; _sigma2A1[ix] = 0.60;
    _A20[ix] = 0.26; _sigma1A2[ix] = 1.32; _sigma2A2[ix] = 0.54;
    _T10[ix] = 0.32; _sigma1T1[ix] = 0.66; _sigma2T1[ix] = 0.31;
    _T20[ix] = 0.32; _sigma1T2[ix] = 0.98; _sigma2T2[ix] = 0.90;
    _T30[ix] = 0.23; _sigma1T3[ix] = 1.42; _sigma2T3[ix] = 0.62;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B_s to eta
  addFormFactor(-531, 221  ,0,-3,5,3);
  _fplus0[39] = 0.36; _sigma1fp[39] = 0.60; _sigma2fp[39] = 0.20;
  _f00[39]    = 0.36; _sigma1f0[39] = 0.80; _sigma2f0[39] = 0.40;
  _fT0[39]    = 0.36; _sigma1fT[39] = 0.58; _sigma2fT[39] = 0.18;
  _massP[39] = 5.37*GeV; _massV[39] = 5.42*GeV;
  // B_s to eta'
  addFormFactor(-531, 331  ,0,-3,5,3);
  _fplus0[40] = 0.36; _sigma1fp[40] = 0.60; _sigma2fp[40] = 0.20;
  _f00[40]    = 0.36; _sigma1f0[40] = 0.80; _sigma2f0[40] = 0.45;
  _fT0[40]    = 0.39; _sigma1fT[40] = 0.58; _sigma2fT[40] = 0.18;
  _massP[40] = 5.37*GeV; _massV[40] = 5.42*GeV;
  // B_s to phi
  addFormFactor(-531, 333  ,1,-3,5,3);
  _V0[41]  = 0.44; _sigma1V0[41] = 0.62; _sigma2V0[41] = 0.20;
  _A00[41] = 0.42; _sigma1A0[41] = 0.55; _sigma2A0[41] = 0.12;
  _A10[41] = 0.34; _sigma1A1[41] = 0.73; _sigma2A1[41] = 0.42;
  _A20[41] = 0.31; _sigma1A2[41] = 1.30; _sigma2A2[41] = 0.52;
  _T10[41] = 0.38; _sigma1T1[41] = 0.62; _sigma2T1[41] = 0.20;
  _T20[41] = 0.38; _sigma1T2[41] = 0.83; _sigma2T2[41] = 0.71;
  _T30[41] = 0.26; _sigma1T3[41] = 1.41; _sigma2T3[41] = 0.57;
  _massP[41] = 5.37*GeV; _massV[41] = 5.42*GeV;
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = 2./9.*Constants::pi;
}

void MelikhovStechFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_fplus0.size()||isize!=_sigma1fp.size()||isize!=_sigma2fp.size()||
     isize!=_f00.size()   ||isize!=_sigma1f0.size()||isize!=_sigma2f0.size()||
     isize!=_fT0.size()   ||isize!=_sigma1fT.size()||isize!=_sigma2fT.size()||
     isize!=_V0.size()    ||isize!=_sigma1V0.size()||isize!=_sigma2V0.size()||
     isize!=_A00.size()   ||isize!=_sigma1A0.size()||isize!=_sigma2A0.size()||
     isize!=_A10.size()   ||isize!=_sigma1A1.size()||isize!=_sigma2A1.size()||
     isize!=_A20.size()   ||isize!=_sigma1A2.size()||isize!=_sigma2A2.size()||
     isize!=_T10.size()   ||isize!=_sigma1T1.size()||isize!=_sigma2T1.size()||
     isize!=_T20.size()   ||isize!=_sigma1T2.size()||isize!=_sigma2T2.size()||
     isize!=_T30.size()   ||isize!=_sigma1T3.size()||isize!=_sigma2T3.size()||
     isize!=_massP.size() ||isize!=_massV.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "MelikhovStechFormFactor::doinit()" 
			  << Exception::abortnow;
}

void MelikhovStechFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fplus0 << _sigma1fp << _sigma2fp << _f00 << _sigma1f0 << _sigma2f0 << _fT0 
     << _sigma1fT << _sigma2fT << _V0 << _sigma1V0 << _sigma2V0 << _A00 << _sigma1A0 
     << _sigma2A0 << _A10 << _sigma1A1 << _sigma2A1 << _A20 << _sigma1A2 << _sigma2A2 
     << _T10 << _sigma1T1 << _sigma2T1 << _T20 << _sigma1T2 << _sigma2T2 << _T30 
     << _sigma1T3 << _sigma2T3 << ounit(_massP,GeV) << ounit(_massV,GeV) << _thetaeta;
}

void MelikhovStechFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fplus0 >> _sigma1fp >> _sigma2fp >> _f00 >> _sigma1f0 >> _sigma2f0 >> _fT0 
     >> _sigma1fT >> _sigma2fT >> _V0 >> _sigma1V0 >> _sigma2V0 >> _A00 >> _sigma1A0 
     >> _sigma2A0 >> _A10 >> _sigma1A1 >> _sigma2A1 >> _A20 >> _sigma1A2 >> _sigma2A2 
     >> _T10 >> _sigma1T1 >> _sigma2T1 >> _T20 >> _sigma1T2 >> _sigma2T2 >> _T30 
     >> _sigma1T3 >> _sigma2T3 >> iunit(_massP,GeV) >> iunit(_massV,GeV) >> _thetaeta;
}

ClassDescription<MelikhovStechFormFactor> MelikhovStechFormFactor::initMelikhovStechFormFactor;
// Definition of the static class description member.

void MelikhovStechFormFactor::Init() {

  static ClassDocumentation<MelikhovStechFormFactor> documentation
    ("The MelikhovStechFormFactor class is the implementation of the"
     " form factors from Phys. Rev. D62  014006 (2000).",
     "The form factors of \\cite{Melikhov:2000yu} were used.",
     "\\bibitem{Melikhov:2000yu} D.~Melikhov and B.~Stech,\n"
     "Phys.\\ Rev.\\  D {\\bf 62} (2000) 014006 [arXiv:hep-ph/0001113].\n"
     "%%CITATION = PHRVA,D62,014006;%%\n");

  static ParVector<MelikhovStechFormFactor,double> interfaceFPlus0
    ("FPlus0",
     "The value of the F_+ form factor at q^2=0",
     &MelikhovStechFormFactor::_fplus0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_1
    ("F+sigma_1",
     "The value of sigma_1 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma1fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_2
    ("F+sigma_2",
     "The value of sigma_2 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma2fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF00
    ("F00",
     "The value of the F_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_f00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_1
    ("F0sigma_1",
     "The value of sigma_1 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma1f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_2
    ("F0sigma_2",
     "The value of sigma_2 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma2f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFT0
    ("FT0",
     "The value of the F_T form factor at q^2=0",
     &MelikhovStechFormFactor::_fT0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_1
    ("FTsigma_1",
     "The value of sigma_1 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma1fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_2
    ("FTsigma_2",
     "The value of sigma_2 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma2fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV00
    ("V00",
     "The value of the V_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_1
    ("V0sigma_1",
     "The value of sigma_1 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma1V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_2
    ("V0sigma_2",
     "The value of sigma_2 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma2V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA00
    ("A00",
     "The value of the A_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_A00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_1
    ("A0sigma_1",
     "The value of sigma_1 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma1A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_2
    ("A0sigma_2",
     "The value of sigma_2 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma2A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA10
    ("A10",
     "The value of the A_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_A10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_1
    ("A1sigma_1",
     "The value of sigma_1 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma1A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_2
    ("A1sigma_2",
     "The value of sigma_2 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma2A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA20
    ("A20",
     "The value of the A_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_A20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_1
    ("A2sigma_1",
     "The value of sigma_1 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma1A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_2
    ("A2sigma_2",
     "The value of sigma_2 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma2A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT10
    ("T10",
     "The value of the T_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_T10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_1
    ("T1sigma_1",
     "The value of sigma_1 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma1T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_2
    ("T1sigma_2",
     "The value of sigma_2 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma2T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT20
    ("T20",
     "The value of the T_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_T20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_1
    ("T2sigma_1",
     "The value of sigma_1 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma1T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_2
    ("T2sigma_2",
     "The value of sigma_2 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma2T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT30
    ("T30",
     "The value of the T_3 form factor at q^2=0",
     &MelikhovStechFormFactor::_T30, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_1
    ("T3sigma_1",
     "The value of sigma_1 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma1T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_2
    ("T3sigma_2",
     "The value of sigma_2 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma2T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassP
    ("MassP",
     "The mass of the pseudoscalar for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massP, GeV, -1, ZERO, ZERO, 10.*GeV,
     false, false, false);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassV
    ("MassV",
     "The mass of the vector for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massV, GeV, -1, ZERO, ZERO, 10.*GeV,
     false, false, false);

  static Parameter< MelikhovStechFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     & MelikhovStechFormFactor::_thetaeta, 2.*Constants::pi/9.,
     -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void MelikhovStechFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
						     int,int id1,
						     Energy, Energy,Complex & f0,
						     Complex & fp) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  fp = _fplus0[mode]/(1.-ratio)/(1.-ratio*(_sigma1fp[mode]-_sigma2fp[mode]*ratio));
  f0 = _f00[mode]              /(1.-ratio*(_sigma1f0[mode]-_sigma2f0[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>3) return;
  double fact;
  if(id1==ParticleID::eta) {
    if(abs(outquark)==3) fact =-sin(_thetaeta);
    else                 fact = cos(_thetaeta)*sqrt(0.5);
  }
  else if(id1==ParticleID::etaprime) {
    if(abs(outquark)==3)  fact=cos(_thetaeta);
    else                  fact=sin(_thetaeta);
  }
  else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
  else                                            fact= sqrt(0.5);
  f0 *= fact;
  fp *= fact;
}
  
void MelikhovStechFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int,int id1,
							  Energy, Energy,
							  Complex & fT) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  fT = _fT0[mode]/(1.-ratio)/(1.-ratio*(_sigma1fT[mode]-_sigma2fT[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>3) return;
  double fact;
  if(id1==ParticleID::eta) {
    if(abs(outquark)==3) fact = -sin(_thetaeta);
    else                 fact =  cos(_thetaeta)*sqrt(0.5);
  }
  else if(id1==ParticleID::etaprime) {
    if(abs(outquark)==3) fact = cos(_thetaeta);
    else                 fact = sin(_thetaeta);
  }
  else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
  else                                            fact= sqrt(0.5);
  fT*=fact;
}

void MelikhovStechFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int mode,
						     int, int id1,
						     Energy, Energy,Complex & A0,
						     Complex & A1,Complex & A2,
						     Complex & V) const {
  useMe();
  double ratio(q2/sqr(_massV[mode])),ratioP(q2/sqr(_massP[mode]));
  A0= _A00[mode]/(1.-ratioP)/(1.-ratioP*(_sigma1A0[mode]-_sigma2A0[mode]*ratioP));
  A1= _A10[mode]            /(1.-ratio *(_sigma1A1[mode]-_sigma2A1[mode]*ratio));
  A2= _A20[mode]            /(1.-ratio *(_sigma1A2[mode]-_sigma2A2[mode]*ratio));
  V =-_V0[mode] /(1.-ratio )/(1.-ratio *(_sigma1V0[mode]-_sigma2V0[mode]*ratio ));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>2) return;
  double fact(sqrt(0.5));
  if(id1==ParticleID::rho0&&abs(outquark)==1) fact=-fact;
  A0 *= fact;
  A1 *= fact;
  A2 *= fact;
  V  *= fact;
}

void MelikhovStechFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int,int id1,
							  Energy, Energy,
							  Complex & T1,Complex & T2,
							  Complex & T3) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  T1= _T10[mode]/(1.-ratio)/(1.-ratio*(_sigma1T1[mode]-_sigma2T1[mode]*ratio));
  T2= _T20[mode]           /(1.-ratio*(_sigma1T2[mode]-_sigma2T2[mode]*ratio));
  T3= _T30[mode]           /(1.-ratio*(_sigma1T3[mode]-_sigma2T3[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>2) return;
  double fact(sqrt(0.5));
  if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
  T1 *= fact;
  T2 *= fact;
  T3 *= fact;
}

void MelikhovStechFormFactor::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::MelikhovStechFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":FPlus0 " 
	     << ix << "  " << _fplus0[ix] << "\n";
      output << "newdef " << name() << ":F+sigma_1 " 
	     << ix << "  " << _sigma1fp[ix] << "\n";
      output << "newdef " << name() << ":F+sigma_2 " 
	     << ix << "  " << _sigma2fp[ix] << "\n";
      output << "newdef " << name() << ":F00 " 
	     << ix << "  " << _f00[ix] << "\n";
      output << "newdef " << name() << ":F0sigma_1 " 
	     << ix << "  " << _sigma1f0[ix] << "\n";
      output << "newdef " << name() << ":F0sigma_2 " 
	     << ix << "  " << _sigma2f0[ix] << "\n";
      output << "newdef " << name() << ":FT0 " 
	     << ix << "  " << _fT0[ix] << "\n";
      output << "newdef " << name() << ":FTsigma_1 " 
	     << ix << "  " << _sigma1fT[ix] << "\n";
      output << "newdef " << name() << ":FTsigma_2 " 
	     << ix << "  " << _sigma2fT[ix] << "\n";
      output << "newdef " << name() << ":V00 " 
	     << ix << "  " << _V0[ix] << "\n";
      output << "newdef " << name() << ":V0sigma_1 " 
	     << ix << "  " << _sigma1V0[ix] << "\n";
      output << "newdef " << name() << ":V0sigma_2 " 
	     << ix << "  " << _sigma2V0[ix] << "\n";
      output << "newdef " << name() << ":A00 " 
	     << ix << "  " << _A00[ix] << "\n";
      output << "newdef " << name() << ":A0sigma_1 " 
	     << ix << "  " << _sigma1A0[ix] << "\n";
      output << "newdef " << name() << ":A0sigma_2 " 
	     << ix << "  " << _sigma2A0[ix] << "\n";
      output << "newdef " << name() << ":A10 " 
	     << ix << "  " << _A10[ix] << "\n";
      output << "newdef " << name() << ":A1sigma_1 " 
	     << ix << "  " << _sigma1A1[ix] << "\n";
      output << "newdef " << name() << ":A1sigma_2 " 
	     << ix << "  " << _sigma2A1[ix] << "\n";
      output << "newdef " << name() << ":A20 " 
	     << ix << "  " << _A20[ix] << "\n";
      output << "newdef " << name() << ":A2sigma_1 " 
	     << ix << "  " << _sigma1A2[ix] << "\n";
      output << "newdef " << name() << ":A2sigma_2 " 
	     << ix << "  " << _sigma2A2[ix] << "\n";
      output << "newdef " << name() << ":T10 " 
	     << ix << "  " << _T10[ix] << "\n";
      output << "newdef " << name() << ":T1sigma_1 " 
	     << ix << "  " << _sigma1T1[ix] << "\n";
      output << "newdef " << name() << ":T1sigma_2 " 
	     << ix << "  " << _sigma2T1[ix] << "\n";
      output << "newdef " << name() << ":T20 " 
	     << ix << "  " << _T20[ix] << "\n";
      output << "newdef " << name() << ":T2sigma_1 " 
	     << ix << "  " << _sigma1T2[ix] << "\n";
      output << "newdef " << name() << ":T2sigma_2 " 
	     << ix << "  " << _sigma2T2[ix] << "\n";
      output << "newdef " << name() << ":T30 " 
	     << ix << "  " << _T30[ix] << "\n";
      output << "newdef " << name() << ":T3sigma_1 " 
	     << ix << "  " << _sigma1T3[ix] << "\n";
      output << "newdef " << name() << ":T3sigma_2 " 
	     << ix << "  " << _sigma2T3[ix] << "\n";
      output << "newdef " << name() << ":MassP " 
	     << ix << "  " << _massP[ix]/GeV << "\n";
      output << "newdef " << name() << ":MassV " 
	     << ix << "  " << _massV[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":FPlus0 " 
	     << ix << "  " << _fplus0[ix] << "\n";
      output << "insert " << name() << ":F+sigma_1 " 
	     << ix << "  " << _sigma1fp[ix] << "\n";
      output << "insert " << name() << ":F+sigma_2 " 
	     << ix << "  " << _sigma2fp[ix] << "\n";
      output << "insert " << name() << ":F00 " 
	     << ix << "  " << _f00[ix] << "\n";
      output << "insert " << name() << ":F0sigma_1 " 
	     << ix << "  " << _sigma1f0[ix] << "\n";
      output << "insert " << name() << ":F0sigma_2 " 
	     << ix << "  " << _sigma2f0[ix] << "\n";
      output << "insert " << name() << ":FT0 " 
	     << ix << "  " << _fT0[ix] << "\n";
      output << "insert " << name() << ":FTsigma_1 " 
	     << ix << "  " << _sigma1fT[ix] << "\n";
      output << "insert " << name() << ":FTsigma_2 " 
	     << ix << "  " << _sigma2fT[ix] << "\n";
      output << "insert " << name() << ":V00 " 
	     << ix << "  " << _V0[ix] << "\n";
      output << "insert " << name() << ":V0sigma_1 " 
	     << ix << "  " << _sigma1V0[ix] << "\n";
      output << "insert " << name() << ":V0sigma_2 " 
	     << ix << "  " << _sigma2V0[ix] << "\n";
      output << "insert " << name() << ":A00 " 
	     << ix << "  " << _A00[ix] << "\n";
      output << "insert " << name() << ":A0sigma_1 " 
	     << ix << "  " << _sigma1A0[ix] << "\n";
      output << "insert " << name() << ":A0sigma_2 " 
	     << ix << "  " << _sigma2A0[ix] << "\n";
      output << "insert " << name() << ":A10 " 
	     << ix << "  " << _A10[ix] << "\n";
      output << "insert " << name() << ":A1sigma_1 " 
	     << ix << "  " << _sigma1A1[ix] << "\n";
      output << "insert " << name() << ":A1sigma_2 " 
	     << ix << "  " << _sigma2A1[ix] << "\n";
      output << "insert " << name() << ":A20 " 
	     << ix << "  " << _A20[ix] << "\n";
      output << "insert " << name() << ":A2sigma_1 " 
	     << ix << "  " << _sigma1A2[ix] << "\n";
      output << "insert " << name() << ":A2sigma_2 " 
	     << ix << "  " << _sigma2A2[ix] << "\n";
      output << "insert " << name() << ":T10 " 
	     << ix << "  " << _T10[ix] << "\n";
      output << "insert " << name() << ":T1sigma_1 " 
	     << ix << "  " << _sigma1T1[ix] << "\n";
      output << "insert " << name() << ":T1sigma_2 " 
	     << ix << "  " << _sigma2T1[ix] << "\n";
      output << "insert " << name() << ":T20 " 
	     << ix << "  " << _T20[ix] << "\n";
      output << "insert " << name() << ":T2sigma_1 " 
	     << ix << "  " << _sigma1T2[ix] << "\n";
      output << "insert " << name() << ":T2sigma_2 " 
	     << ix << "  " << _sigma2T2[ix] << "\n";
      output << "insert " << name() << ":T30 " 
	     << ix << "  " << _T30[ix] << "\n";
      output << "insert " << name() << ":T3sigma_1 " 
	     << ix << "  " << _sigma1T3[ix] << "\n";
      output << "insert " << name() << ":T3sigma_2 " 
	     << ix << "  " << _sigma2T3[ix] << "\n";
      output << "insert " << name() << ":MassP " 
	     << ix << "  " << _massP[ix]/GeV << "\n";
      output << "insert " << name() << ":MassV " 
	     << ix << "  " << _massV[ix]/GeV << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
