// -*- C++ -*-
//
// HiddenValleyAlpha.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenValleyAlpha class.
//
#include "HiddenValleyAlpha.h"
#include "HiddenValleyModel.h"
#include "EnumParticles.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

IBPtr HiddenValleyAlpha::clone() const {
  return new_ptr(*this);
}

IBPtr HiddenValleyAlpha::fullclone() const {
  return new_ptr(*this);
}

void HiddenValleyAlpha::persistentOutput(PersistentOStream & os) const {
  os << _asType << _asMaxNP << ounit(_qmin,GeV) << _nloop << _thresopt 
     << ounit(_lambdain,GeV) << _alphain 
     << _tolerance << _maxtry << _alphamin
     << ounit(_thresholds,GeV) << ounit(_lambda,GeV) << _ca << _cf << _tr;
}

void HiddenValleyAlpha::persistentInput(PersistentIStream & is, int) {
  is >> _asType >> _asMaxNP >> iunit(_qmin,GeV) >> _nloop >> _thresopt
     >> iunit(_lambdain,GeV) >> _alphain 
     >> _tolerance >> _maxtry >> _alphamin
     >> iunit(_thresholds,GeV) >> iunit(_lambda,GeV) >> _ca >> _cf >> _tr;
}

ClassDescription<HiddenValleyAlpha> HiddenValleyAlpha::initHiddenValleyAlpha;
// Definition of the static class description member.

void HiddenValleyAlpha::Init() {

  static ClassDocumentation<HiddenValleyAlpha> documentation
    ("This (concrete) class describes the QCD alpha running.");

  static Switch<HiddenValleyAlpha, int> intAsType
    ("NPAlphaS",
     "Behaviour of AlphaS in the NP region",
     &HiddenValleyAlpha::_asType, 1, false, false);
  static SwitchOption intAsTypeZero
    (intAsType, "Zero","zero below Q_min", 1);
  static SwitchOption intAsTypeConst
    (intAsType, "Const","const as(qmin) below Q_min", 2);
  static SwitchOption intAsTypeLin
    (intAsType, "Linear","growing linearly below Q_min", 3);
  static SwitchOption intAsTypeQuad
    (intAsType, "Quadratic","growing quadratically below Q_min", 4);
  static SwitchOption intAsTypeExx1
    (intAsType, "Exx1", "quadratic from AlphaMaxNP down to as(Q_min)", 5);
  static SwitchOption intAsTypeExx2
    (intAsType, "Exx2", "const = AlphaMaxNP below Q_min", 6);

  // default such that as(qmin) = 1 in the current parametrization.
  // min = Lambda3
  static Parameter<HiddenValleyAlpha,Energy> intQmin
    ("Qmin", "Q < Qmin is treated with NP parametrization as of (unit [GeV])",
     &HiddenValleyAlpha::_qmin, GeV, 0.630882*GeV, 0.330445*GeV,
     100.0*GeV,false,false,false);

  static Parameter<HiddenValleyAlpha,double> interfaceAlphaMaxNP
    ("AlphaMaxNP",
     "Max value of alpha in NP region, only relevant if NPAlphaS = 5,6",
     &HiddenValleyAlpha::_asMaxNP, 1.0, 0., 100.0,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,unsigned int> interfaceNumberOfLoops
    ("NumberOfLoops",
     "The number of loops to use in the alpha_S calculation",
     &HiddenValleyAlpha::_nloop, 2, 1, 2,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "Input value of Lambda_MSBar",
     &HiddenValleyAlpha::_lambdain, MeV, 0.208364*GeV, 100.0*MeV, 500.0*MeV,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,double> interfaceAlphaMZ
    ("AlphaMZ",
     "The input value of the strong coupling at the Z mass ",
     &HiddenValleyAlpha::_alphain, 0.118, 0.1, 0.2,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,double> interfaceTolerance
    ("Tolerance",
     "The tolerance for discontinuities in alphaS at thresholds.",
     &HiddenValleyAlpha::_tolerance, 1e-10, 1e-20, 1e-4,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,unsigned int> interfaceMaximumIterations
    ("MaximumIterations",
     "The maximum number of iterations for the Newton-Raphson method to converge.",
     &HiddenValleyAlpha::_maxtry, 100, 10, 1000,
     false, false, Interface::limited);

  static Switch<HiddenValleyAlpha,bool> interfaceThresholdOption
    ("ThresholdOption",
     "Whether to use the consistuent or normal masses for the thresholds",
     &HiddenValleyAlpha::_thresopt, true, false, false);
  static SwitchOption interfaceThresholdOptionCurrent
    (interfaceThresholdOption,
     "Current",
     "Use the current masses",
     true);
  static SwitchOption interfaceThresholdOptionConstituent
    (interfaceThresholdOption,
     "Constituent",
     "Use the constitent masses.",
     false);

}

double HiddenValleyAlpha::value(const Energy2 scale) const {
  pair<short,Energy> nflam;
  Energy q = sqrt(scale);
  double val(0.);
  // special handling if the scale is less than Qmin
  if (q < _qmin) {
    nflam = getLamNfTwoLoop(_qmin); 
    double val0 = alphaS(_qmin, nflam.second, nflam.first);
    switch (_asType) {
    case 1: 
      // flat, zero; the default type with no NP effects.
      val = 0.; 
      break; 
    case 2: 
      // flat, non-zero alpha_s = alpha_s(q2min).
      val = val0;
      break; 
    case 3: 
      // linear in q
      val = val0*q/_qmin;
      break; 
    case 4:
      // quadratic in q
      val = val0*sqr(q/_qmin);
      break; 
    case 5:
      // quadratic in q, starting off at asMaxNP, ending on as(qmin)
      val = (val0 - _asMaxNP)*sqr(q/_qmin) + _asMaxNP;
      break; 
    case 6:
      // just asMaxNP and constant
      val = _asMaxNP;
      break; 
    }
  } 
  else {
    // the 'ordinary' case    
    nflam = getLamNfTwoLoop(q); 
    val = alphaS(q, nflam.second, nflam.first);
  }
  return scaleFactor() * val;
}

double HiddenValleyAlpha::overestimateValue() const {
  return scaleFactor() * _alphamin; 
}

double HiddenValleyAlpha::ratio(const Energy2 scale) const {
  pair<short,Energy> nflam;
  Energy q = sqrt(scale);
  double val(0.);
  // special handling if the scale is less than Qmin
  if (q < _qmin) {
    nflam = getLamNfTwoLoop(_qmin); 
    double val0 = alphaS(_qmin, nflam.second, nflam.first);
    switch (_asType) {
    case 1: 
      // flat, zero; the default type with no NP effects.
      val = 0.; 
      break; 
    case 2: 
      // flat, non-zero alpha_s = alpha_s(q2min).
      val = val0;
      break; 
    case 3: 
      // linear in q
      val = val0*q/_qmin;
      break; 
    case 4:
      // quadratic in q
      val = val0*sqr(q/_qmin);
      break; 
    case 5:
      // quadratic in q, starting off at asMaxNP, ending on as(qmin)
      val = (val0 - _asMaxNP)*sqr(q/_qmin) + _asMaxNP;
      break; 
    case 6:
      // just asMaxNP and constant
      val = _asMaxNP;
      break; 
    }
  } else {
    // the 'ordinary' case    
    nflam = getLamNfTwoLoop(q); 
    val = alphaS(q, nflam.second, nflam.first);
  }
  // denominator 
  return val/_alphamin;  
}

Energy HiddenValleyAlpha::computeLambda(Energy match,
				     double alpha,
				     unsigned int nflav) const {
  Energy lamtest=200.0*MeV;
  double xtest;
  unsigned int ntry=0;
  do {
    ++ntry;
    xtest=log(sqr(match/lamtest));
    xtest+= (alpha-alphaS(match,lamtest,nflav))/derivativealphaS(match,lamtest,nflav);
    lamtest=match/exp(0.5*xtest);
  }
  while(abs(alpha-alphaS(match,lamtest,nflav)) > _tolerance && ntry < _maxtry);
  return lamtest;
}

pair<short, Energy> HiddenValleyAlpha::getLamNfTwoLoop(Energy q) const {
  unsigned int ix=1;
  for(;ix<_thresholds.size();++ix) {
    if(q<_thresholds[ix]) break;
  }
  if(ix==_thresholds.size()) --ix;
  --ix;
  return pair<short,Energy>(ix, _lambda[ix]);
}

void HiddenValleyAlpha::doinit() {
  ShowerAlpha::doinit();
  // get the model for parameters
  tcHiddenValleyPtr model = dynamic_ptr_cast<tcHiddenValleyPtr>
    (generator()->standardModel());
  if(!model) throw InitException() << "Must be using the HiddenValleyModel"
				   << " in HiddenValleyAlpha::doinit()" 
				   << Exception::runerror;
  // get the colour factors
  _ca = model->CA();
  _cf = model->CF();
  _tr = model->TR();
  // get the thresholds
  _thresholds.push_back(_qmin);
  for(unsigned int ix=1;ix<=model->NF();++ix) {
    _thresholds.push_back(getParticleData(HiddenID::darkGluon+long(ix))->mass());
  }
  _lambda.resize(_thresholds.size());
  Energy mz = getParticleData(ThePEG::ParticleID::Z0)->mass();
  unsigned int nf;
  for(nf=0;nf<_thresholds.size();++nf) {
    if(mz<_thresholds[nf]) break;
  }
  nf-=1;
  // value of Lambda from alphas if needed using Newton-Raphson
  _lambda[nf] = computeLambda(mz,_alphain,nf-1);
  // compute the threshold matching
  // above the Z mass
  for(int ix=nf;ix<int(_thresholds.size())-1;++ix) {
    _lambda[ix+1] = computeLambda(_thresholds[ix+1],alphaS(_thresholds[ix+1],
							   _lambda[ix],ix),ix+1);
  }
  // below Z mass
  for(int ix=nf-1;ix>=0;--ix) {
    _lambda[ix] = computeLambda(_thresholds[ix+1],alphaS(_thresholds[ix+1],
							 _lambda[ix+1],ix+1),ix);
  }
  // compute the maximum value of as 
  if ( _asType < 5 ) _alphamin = value(sqr(_qmin)+1.0e-8*sqr(MeV)); // approx as = 1
  else _alphamin = max(_asMaxNP, value(sqr(_qmin)+1.0e-8*sqr(MeV))); 
  // check consistency lambda_3 < qmin
  if(_lambda[0]>_qmin)
    Throw<InitException>() << "The value of Qmin is less than Lambda in "
			   << _qmin/GeV << " < " << _lambda[0]/GeV 
			   << " HiddenValleyAlpha::doinit " << Exception::abortnow;
  for(Energy scale=_qmin;scale<200.*GeV;scale+=0.2*GeV)
    cerr << scale/GeV << "\t" << value(sqr(scale)) << "\n";
}

double HiddenValleyAlpha::derivativealphaS(Energy q, Energy lam, int nf) const {
  using Constants::pi;
  double lx = log(sqr(q/lam));
  // N.B. b_1 is divided by 2 due Hw++ convention
  double b0 = 11./3.*_ca - 4./3.*_tr*nf;
  double b1 = 17./3.*sqr(_ca) - nf*_tr*(10./3.*_ca+2.*_cf);
  if(_nloop==1)
    return -4.*pi/(b0*sqr(lx));
  else if(_nloop==2)
    return -4.*pi/(b0*sqr(lx))*(1.-2.*b1/sqr(b0)/lx*(1.-2.*log(lx)));
  else
    assert(false);
}

double HiddenValleyAlpha::alphaS(Energy q, Energy lam, int nf) const {
  using Constants::pi;
  double lx(log(sqr(q/lam)));
  // N.B. b_1 is divided by 2 due Hw++ convention
  double b0 = 11./3.*_ca - 4./3.*_tr*nf;
  double b1 = 17./3.*sqr(_ca) - nf*_tr*(10./3.*_ca+2.*_cf);
  // one loop
  if(_nloop==1)
    return 4.*pi/(b0*lx);
  // two loop
  else if(_nloop==2) {
    return 4.*pi/(b0*lx)*(1.-2.*b1/sqr(b0)*log(lx)/lx);
  }
  else
    assert(false);
}
