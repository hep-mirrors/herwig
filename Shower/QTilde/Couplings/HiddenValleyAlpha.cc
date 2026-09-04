// -*- C++ -*-
//
// HiddenValleyAlpha.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiddenValleyAlpha class.
//
#include "HiddenValleyAlpha.h"
#include "Herwig/Models/HiddenValley/HiddenValleyModel.h"
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
     << ounit(_lambdain,GeV) << _tolerance << _maxtry << _alphamin
     << ounit(_thresholds,GeV) << ounit(_lambda,GeV) << _ca << _cf << _tr << _nf_light;
}

void HiddenValleyAlpha::persistentInput(PersistentIStream & is, int) {
  is >> _asType >> _asMaxNP >> iunit(_qmin,GeV) >> _nloop >> _thresopt
     >> iunit(_lambdain,GeV) >> _tolerance >> _maxtry >> _alphamin
     >> iunit(_thresholds,GeV) >> iunit(_lambda,GeV) >> _ca >> _cf >> _tr >> _nf_light;
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

  static Switch<HiddenValleyAlpha,bool> interfaceLambdaOption
    ("LambdaOption",
     "Option for the calculation of the Lambda used in the simulation from the input"
     " Lambda_MSbar",
     &HiddenValleyAlpha::_lambdaopt, false, false, false);
  static SwitchOption interfaceLambdaOptionfalse
    (interfaceLambdaOption,
     "Same",
     "Use the same value",
     false);
  static SwitchOption interfaceLambdaOptionConvert
    (interfaceLambdaOption,
     "Convert",
     "Use the conversion to the Herwig scheme from NPB349, 635",
     true);

  static Parameter<HiddenValleyAlpha,Energy> interfaceLambdaDark
    ("LambdaDark",
     "Input value of Lambda_MSBar",
     &HiddenValleyAlpha::_lambdain, GeV, 0.208364*GeV, 100.0*MeV, 100.0*GeV,
     false, false, Interface::limited);

  static Parameter<HiddenValleyAlpha,double> interfaceAlphaIn
    ("AlphaIn",
     "The input value of the dark coupling at the chosen InputScale (default: MZ)",
     &HiddenValleyAlpha::_alphain, 0.118, 0.1, 0.2,
     false, false, Interface::limited);
  
  static Parameter<HiddenValleyAlpha,Energy> interfaceInputScale
    ("InputScale",
     "An optional input scale. MZ will be used if not set.",
     &HiddenValleyAlpha::_optInputScale, GeV, 91.1876_GeV, ZERO, 0*GeV,
     false, false, Interface::lowerlim);

static Switch<HiddenValleyAlpha,bool> interfaceInputOption
    ("InputOption",
     "Option for inputing the initial value of the coupling",
     &HiddenValleyAlpha::_inopt, false, false, false);
  static SwitchOption interfaceInputOptionAlphaMZ
    (interfaceInputOption,
     "AlphaIn",
     "Use the value of alpha at the InputScale to calculate the coupling",
     true);
  static SwitchOption interfaceInputOptionLambdaQCD
    (interfaceInputOption,
     "LambdaDark",
     "Use the input value of Lambda to calculate the coupling",
     false);
  
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
  
  static Command<HiddenValleyAlpha> interfacecheck
    ("check",
     "check",
     &HiddenValleyAlpha::check, false);
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

double HiddenValleyAlpha::ratio(const Energy2 scale,double) const {
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
  for(;ix<_thresholds.size();ix++) {
    if(q<_thresholds[ix]) break;
  }
  --ix;
  return pair<short,Energy>(ix+_nf_light, _lambda[ix]);
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
  _nf_light = 0;
  for(unsigned int ix=1;ix<=model->NF();++ix) {
    Energy qmass;
    if(_thresopt)
      qmass = getParticleData(ParticleID::darkg+long(ix))->mass();
    else
      qmass = getParticleData(ParticleID::darkg+long(ix))->constituentMass();
    if (qmass > _qmin) _thresholds.push_back(qmass);
    else _nf_light++;
  }
  _lambda.resize(_thresholds.size());


  unsigned int nf = _nf_light;

  
  if(_inopt) {
    Energy refscale = _optInputScale == ZERO ? getParticleData(ParticleID::Z0)->mass() : _optInputScale;
    if (refscale < _qmin) throw InitException() << "InputScale for Alpha_Dark must be higher than Qmin"
				   << Exception::runerror;
    unsigned int ix;
    for(ix=0;ix<_thresholds.size();++ix) {
        if(refscale<_thresholds[ix]) break;
    }
    ix-=1;
    nf += ix;
    _lambda[ix]=computeLambda(refscale,
			     _alphain,nf);
  }
  // otherwise it was an input parameter, take this as the value below heavy quark thresholds
  else{_lambda[0]=_lambdain;}
  // convert lambda to the Monte Carlo scheme if needed
  using Constants::pi;
  if(_lambdaopt) _lambda[0] *=exp((0.5*_ca*(67./3.-sqr(pi))-_tr*nf*10./3.)/(11*_ca-2*nf))/sqrt(2.);

  // Compute the threshold matching
  // above the reference scale
  for(int ix=nf-_nf_light;ix<int(_thresholds.size())-1;++ix) {
    _lambda[ix+1] = computeLambda(_thresholds[ix+1],alphaS(_thresholds[ix+1],
							   _lambda[ix],ix+_nf_light),ix+_nf_light+1);
  }
  // below the reference scale
  for(int ix=nf-1-_nf_light;ix>=0;--ix) {
    _lambda[ix] = computeLambda(_thresholds[ix+1],alphaS(_thresholds[ix+1],
							 _lambda[ix+1],ix+_nf_light+1),ix+_nf_light);
  }
  // compute the maximum value of as
  if ( _asType < 5 ) _alphamin = value(sqr(_qmin)); // approx as = 1
  else _alphamin = max(_asMaxNP, value(sqr(_qmin)));
  // check consistency lambda_3 < qmin
  if(_lambda[0]>_qmin)
    Throw<InitException>() << "The value of Qmin is less than Lambda in "
			   << _qmin/GeV << " < " << _lambda[0]/GeV
			   << " HiddenValleyAlpha::doinit " << Exception::abortnow;
  if (_do_check) {
      
    Repository::clog() << "checking alpha_dark in range [" << _check_Q_low << "," << _check_Q_high << "] GeV in "
	            	   << _check_n_steps << " steps.\nResults are written to " << _check_fname << "\n";

    double step_width = (_check_Q_high-_check_Q_low)/_check_n_steps;

    ofstream out (_check_fname.c_str());

    for (long k = 0; k <= _check_n_steps; ++k) {

      Energy Q = _check_Q_low*GeV + k*step_width*GeV;

      out << (Q/GeV) << " " << value(Q*Q) << "\n";

    }
  }
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
    return -4.*pi/(b0*sqr(lx))*(1.+2.*b1/sqr(b0)/lx*(1.-2.*log(lx)));
  else
    assert(false);
  return 0.;
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
  return 0.;
}

string HiddenValleyAlpha::check(string args) {

  istringstream argin(args);

  _do_check = true;

  argin >> _check_Q_low >> _check_Q_high >> _check_n_steps;

  argin >> _check_fname;
  return "Check flags set";
}
