// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FortranAlphaQCD class.
//
#include "FortranAlphaQCD.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

void FortranAlphaQCD::persistentOutput(PersistentOStream & os) const {
  //os << _thresholds << _lambdain << _lambda5 << _lambda3 << _maxtry 
  //   << _alphamax << _match << _ccoeff << _bcoeff;
}

void FortranAlphaQCD::persistentInput(PersistentIStream & is, int) {
  //is >> _thresholds >> _lambdain >> _lambda5 >> _lambda3 >> _maxtry 
  //   >> _alphamax >> _match >> _ccoeff >> _bcoeff ;
}

ClassDescription<FortranAlphaQCD> FortranAlphaQCD::initFortranAlphaQCD;
// Definition of the static class description member.

void FortranAlphaQCD::Init() {

  static ClassDocumentation<FortranAlphaQCD> documentation
    ("This (concrete) class describes the QCD alpha running.");

  static Parameter<FortranAlphaQCD,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "Input value of Lambda_MSBar",
     &FortranAlphaQCD::_lambdain, MeV, 0.180*GeV, 100.0*MeV, 500.0*MeV,
     false, false, Interface::limited);

  static Parameter<FortranAlphaQCD,unsigned int> interfaceMaximumIterations
    ("MaximumIterations",
     "The maximum number of iterations for the Newton-Raphson method to converge.",
     &FortranAlphaQCD::_maxtry, 100, 10, 1000,
     false, false, Interface::limited);

  static Parameter<FortranAlphaQCD,double> interfaceMaximumValue
    ("MaximumValue",
     "The maximum value of alpha_S",
     &FortranAlphaQCD::_alphamax, 1.0, 0.1, 10.0,
     false, false, Interface::limited);
}

void FortranAlphaQCD::doinit() throw(InitException) {
  ShowerAlpha::doinit();
  _lambda5=_lambdain*exp(0.5*(67.-3.*sqr(pi)-50./3.)/23.)/sqrt(2.);
  // thresholds
  _thresholds.resize(4);
  for(int ix=1;ix<4;++ix) _thresholds[ix]=getParticleData(ix+3)->constituentMass();
  // beta function coefficients
  double ca=3.,cf=4./3.;
  _bcoeff.resize(4);
  _ccoeff.resize(4);
  for(unsigned int ix=0;ix<4;++ix) {
    _bcoeff[ix]=(11.*ca- 2.*(ix+3))/(12.*pi);
    _ccoeff[ix]=(17.*sqr(ca)-(5.*ca+3.*cf)*(ix+3))/(24.*sqr(pi))/sqr(_bcoeff[ix]);
  }
  // calculate threshold matching
  double rho=2.*log(_thresholds[3]/_lambda5);
  double rat=log(rho)/rho;
  _match.resize(3);
  _match[2]=(_bcoeff[2]/(1.-_ccoeff[2]*rat)-_bcoeff[4]/(1.-_ccoeff[3]*rat))*rho;
  rho=2.*log(_thresholds[2]/_lambda5);
  rat=log(rho)/rho;
  _match[1]=(_bcoeff[2]/(1.-_ccoeff[2]*rat)-_bcoeff[1]/(1.-_ccoeff[1]*rat))*rho;
  rho=2.*log(_thresholds[1]/_lambda5);
  rat=log(rho)/rho;
  _match[0]=(_bcoeff[1]/(1.-_ccoeff[1]*rat)-_bcoeff[0]/(1.-_ccoeff[0]*rat))*rho+_match[1];
  // compute Lambda_3
  double d35=-1./(_bcoeff[0]*_match[0]),drh,rlf,eps(1e-6);
  unsigned int ix=0;
  do {
    rat=log(d35)/d35;
    rlf=_bcoeff[0]*d35/(1.-_ccoeff[0]*rat);
    drh=_bcoeff[0]*(rlf+_match[0])*sqr(d35)/((1.-2.*_ccoeff[0]*rat+_ccoeff[0]/d35)*sqr(rlf));
    d35=d35-drh;
  }
  while(ix<_maxtry&&abs(drh)>eps*d35);
  _lambda3=_lambda5*exp(0.5*d35);
//   ofstream output("alpha.top");
//   output << "set limits y 0. 1.\n";
//   for(Energy scale=0.5*GeV;scale<100.*GeV;scale+=0.5*GeV) 
//     output << scale/GeV << " " << value(sqr(scale)) << "\n";
//   output << "join red\n";
//   output << "new frame \n";
//   output << "set limits y 0. 1.\n";
//   for(Energy scale=0.5*GeV;scale<100.*GeV;scale+=0.5*GeV) 
//     output << scale/GeV << " " << oneLoopValue(sqr(scale)) << "\n";
//   output << "join red\n";
//   output << "new frame \n";
//   output << "set limits y 0. 1.\n";
//   for(Energy scale=0.5*GeV;scale<100.*GeV;scale+=0.5*GeV) 
//     output << scale/GeV << " " << twoLoopRatio(sqr(scale)) << "\n";
//   output << "join red\n";
//   output.close();
}

double FortranAlphaQCD::value(const Energy2 scale) const {
  Energy q = sqrt(scale);
  double rho=2.*log(q/_lambda5);
  double rat=log(rho)/rho,rlf;
  if(q>_thresholds[3]) {
    rlf=_bcoeff[4]*rho/(1.-_ccoeff[3]*rat)+_match[2];
  }
  else if(q>_thresholds[2]) {
    rlf=_bcoeff[2]*rho/(1.-_ccoeff[2]*rat);
  }
  else if(q>_thresholds[1]) {
    rlf=_bcoeff[1]*rho/(1.-_ccoeff[1]*rat)+_match[1];
  }
  else {
    rlf=_bcoeff[0]*rho/(1.-_ccoeff[0]*rat)+_match[0];
  }
  rlf = max(1./_alphamax,rlf);
  return scaleFactor()/rlf;
}

double FortranAlphaQCD::overestimateValue() const {
  return scaleFactor() * _alphamax; 
}

double FortranAlphaQCD::ratio(const Energy2 scale) const {
  return value(scale)/_alphamax;
}

double FortranAlphaQCD::oneLoopValue(const Energy2 scale) const {
  Energy q=sqrt(scale);
  double rho=2.*log(q/_lambda3);
  if(rho<0.) cerr << "problem in one loop" << endl;
  return scaleFactor()/(_bcoeff[2]*rho);
}

double FortranAlphaQCD::twoLoopRatio(const Energy2 scale) const {
  Energy q = sqrt(scale);
  double rho=2.*log(q/_lambda5);
  double rat=log(rho)/rho,rlf;
  if(q>_thresholds[3]) {
    rlf=_bcoeff[4]*rho/(1.-_ccoeff[3]*rat)+_match[2];
  }
  else if(q>_thresholds[2]) {
    rlf=_bcoeff[2]*rho/(1.-_ccoeff[2]*rat);
  }
  else if(q>_thresholds[1]) {
    rlf=_bcoeff[1]*rho/(1.-_ccoeff[1]*rat)+_match[1];
  }
  else {
    rlf=_bcoeff[0]*rho/(1.-_ccoeff[0]*rat)+_match[0];
  }
  rlf = max(1./_alphamax,rlf);
  return _bcoeff[2]*2.*log(q/_lambda3)/rlf;
}
