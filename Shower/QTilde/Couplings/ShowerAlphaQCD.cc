// -*- C++ -*-
//
// ShowerAlphaQCD.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQCD class.
//
#include "ShowerAlphaQCD.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

DescribeClass<ShowerAlphaQCD,ShowerAlpha>
describeShowerAlphaQCD("Herwig::ShowerAlphaQCD","HwShower.so");

IBPtr ShowerAlphaQCD::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerAlphaQCD::fullclone() const {
  return new_ptr(*this);
}

void ShowerAlphaQCD::persistentOutput(PersistentOStream & os) const {
  os << _asType << _asMaxNP << ounit(_qmin,GeV) << _nloop << _lambdaopt << _thresopt 
     << ounit(_lambdain,GeV) << _alphain << _inopt
     << _tolerance << _maxtry << _alphamin
     << ounit(_thresholds,GeV) << ounit(_lambda,GeV)
     << _val0 << ounit(_optInputScale,GeV) << ounit(_quarkMasses,GeV);
}

void ShowerAlphaQCD::persistentInput(PersistentIStream & is, int) {
  is >> _asType >> _asMaxNP >> iunit(_qmin,GeV) >> _nloop >> _lambdaopt >> _thresopt
     >> iunit(_lambdain,GeV) >> _alphain >> _inopt
     >> _tolerance >> _maxtry >> _alphamin
     >> iunit(_thresholds,GeV) >> iunit(_lambda,GeV)
     >> _val0 >> iunit(_optInputScale,GeV) >> iunit(_quarkMasses,GeV);
}

void ShowerAlphaQCD::Init() {

  static ClassDocumentation<ShowerAlphaQCD> documentation
    ("This (concrete) class describes the QCD alpha running.");

  static Switch<ShowerAlphaQCD, int> intAsType
    ("NPAlphaS",
     "Behaviour of AlphaS in the NP region",
     &ShowerAlphaQCD::_asType, 1, false, false);
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
  static Parameter<ShowerAlphaQCD,Energy> intQmin
    ("Qmin", "Q < Qmin is treated with NP parametrization as of (unit [GeV]),"
     " if negative it is solved for at initialisation such that alpha_S(Qmin)=AlphaMaxNP",
     &ShowerAlphaQCD::_qmin, GeV, 0.630882*GeV, 0.330445*GeV,
     100.0*GeV,false,false,false);

  static Parameter<ShowerAlphaQCD,double> interfaceAlphaMaxNP
    ("AlphaMaxNP",
     "Max value of alpha in NP region, only relevant if NPAlphaS = 5,6",
     &ShowerAlphaQCD::_asMaxNP, 1.0, 0., 100.0,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCD,unsigned int> interfaceNumberOfLoops
    ("NumberOfLoops",
     "The number of loops to use in the alpha_S calculation",
     &ShowerAlphaQCD::_nloop, 3, 1, 3,
     false, false, Interface::limited);

  static Switch<ShowerAlphaQCD,bool> interfaceLambdaOption
    ("LambdaOption",
     "Option for the calculation of the Lambda used in the simulation from the input"
     " Lambda_MSbar",
     &ShowerAlphaQCD::_lambdaopt, false, false, false);
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

  static Parameter<ShowerAlphaQCD,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "Input value of Lambda_MSBar",
     &ShowerAlphaQCD::_lambdain, MeV, 0.208364*GeV, 100.0*MeV, 500.0*MeV,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCD,double> interfaceAlphaMZ
    ("AlphaMZ",
     "The input value of the strong coupling at the Z mass ",
     &ShowerAlphaQCD::_alphain, 0.118, 0.1, 0.2,
     false, false, Interface::limited);

  static Switch<ShowerAlphaQCD,bool> interfaceInputOption
    ("InputOption",
     "Option for inputing the initial value of the coupling",
     &ShowerAlphaQCD::_inopt, true, false, false);
  static SwitchOption interfaceInputOptionAlphaMZ
    (interfaceInputOption,
     "AlphaMZ",
     "Use the value of alpha at MZ to calculate the coupling",
     true);
  static SwitchOption interfaceInputOptionLambdaQCD
    (interfaceInputOption,
     "LambdaQCD",
     "Use the input value of Lambda to calculate the coupling",
     false);

  static Parameter<ShowerAlphaQCD,double> interfaceTolerance
    ("Tolerance",
     "The tolerance for discontinuities in alphaS at thresholds.",
     &ShowerAlphaQCD::_tolerance, 1e-10, 1e-20, 1e-4,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCD,unsigned int> interfaceMaximumIterations
    ("MaximumIterations",
     "The maximum number of iterations for the Newton-Raphson method to converge.",
     &ShowerAlphaQCD::_maxtry, 100, 10, 1000,
     false, false, Interface::limited);

  static Switch<ShowerAlphaQCD,bool> interfaceThresholdOption
    ("ThresholdOption",
     "Whether to use the consistuent or normal masses for the thresholds",
     &ShowerAlphaQCD::_thresopt, true, false, false);
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

  static Command<ShowerAlphaQCD> interfaceValue
    ("Value",
     "",
     &ShowerAlphaQCD::value, false);

  static Command<ShowerAlphaQCD> interfacecheck
    ("check",
     "check",
     &ShowerAlphaQCD::check, false);

  static Parameter<ShowerAlphaQCD,Energy> interfaceInputScale
    ("InputScale",
     "An optional input scale. MZ will be used if not set.",
     &ShowerAlphaQCD::_optInputScale, GeV, ZERO, ZERO, 0*GeV,
     false, false, Interface::lowerlim);

  static ParVector<ShowerAlphaQCD,Energy> interfaceQuarkMasses
    ("QuarkMasses",
     "The quark masses to be used instead of the masses set in the particle data.",
     &ShowerAlphaQCD::_quarkMasses, GeV, -1, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

void ShowerAlphaQCD::doinit() {
  ShowerAlpha::doinit();
  // calculate the value of 5-flavour lambda 
  // evaluate the initial
  // value of Lambda from alphas if needed using Newton-Raphson
  if(_inopt) {
    _lambda[2]=computeLambda(_optInputScale == ZERO ? 
			     getParticleData(ParticleID::Z0)->mass() :
			     _optInputScale,
			     _alphain,5);
  }
  // otherwise it was an input parameter
  else{_lambda[2]=_lambdain;}
  // convert lambda to the Monte Carlo scheme if needed
  using Constants::pi;
  if(_lambdaopt) _lambda[2] *=exp(0.5*(67.-3.*sqr(pi)-50./3.)/23.)/sqrt(2.);
  // compute the threshold matching
  // top threshold
  for(int ix=1;ix<4;++ix) {
    if ( _quarkMasses.empty() ) {
      if(_thresopt)
	_thresholds[ix]=getParticleData(ix+3)->mass();
      else
	_thresholds[ix]=getParticleData(ix+3)->constituentMass();
    } else {
      // starting at zero rather than one, cf the other alphas's
      _thresholds[ix] = _quarkMasses[ix+2];
    }
  }
  // compute 6 flavour lambda by matching at top mass using Newton Raphson
  _lambda[3]=computeLambda(_thresholds[3],alphaS(_thresholds[3],_lambda[2],5),6);
  // bottom threshold
  // compute 4 flavour lambda by matching at bottom mass using Newton Raphson
  _lambda[1]=computeLambda(_thresholds[2],alphaS(_thresholds[2],_lambda[2],5),4);
  // charm threshold
  // compute 3 flavour lambda by matching at charm mass using Newton Raphson
  _lambda[0]=computeLambda(_thresholds[1],alphaS(_thresholds[1],_lambda[1],4),3);
  // if qmin less than zero solve for alphaS(_qmin) = 1.
  if(_qmin<ZERO ) {
    Energy low  = _lambda[0]+MeV;
    if(value(sqr(low))<_asMaxNP)
      Throw<InitException>() << "The value of Qmin is less than Lambda_3 in"
			     << " ShowerAlphaQCD::doinit " << Exception::abortnow;
    Energy high = 10*GeV;
    while (value(sqr(high))>_asMaxNP) high *=2.;
    double as;
    do {
      _qmin = 0.5*(low+high);
      as  = value(sqr(_qmin));
      if(_asMaxNP>as) high = _qmin;
      else if(_asMaxNP<as) low = _qmin;
    }
    while ( abs(as-_asMaxNP) > _tolerance );
  }
  // final threshold is qmin
  _thresholds[0]=_qmin;
  // value of alphaS at threshold
  pair<short,Energy> nflam = getLamNfTwoLoop(_qmin);
  _val0 = alphaS(_qmin, nflam.second, nflam.first);
  // compute the maximum value of as 
  if ( _asType < 5 )
    _alphamin = _val0;
  else
    _alphamin = max(_asMaxNP, _val0);
  // check consistency lambda_3 < qmin
  if(_lambda[0]>_qmin)
    Throw<InitException>() << "The value of Qmin is less than Lambda_3 in"
			   << " ShowerAlphaQCD::doinit " << Exception::abortnow;
}

string ShowerAlphaQCD::check(string args) {

  doinit();

  istringstream argin(args);

  double Q_low, Q_high;
  long n_steps;

  argin >> Q_low >> Q_high >> n_steps;

  string fname;
  argin >> fname;

  Repository::clog() << "checking alpha_s in range [" << Q_low << "," << Q_high << "] GeV in "
		     << n_steps << " steps.\nResults are written to " << fname << "\n";

  double step_width = (Q_high-Q_low)/n_steps;

  ofstream out (fname.c_str());

  for (long k = 0; k <= n_steps; ++k) {

    Energy Q = Q_low*GeV + k*step_width*GeV;

    out << (Q/GeV) << " " << value(Q*Q) << "\n";

  }

  return "alpha_s check finished";

}


double ShowerAlphaQCD::value(const Energy2 scale) const {
  Energy q = scaleFactor()*sqrt(scale);
  double val(0.);
  // normal case
  if (q >= _qmin) {
    pair<short,Energy> nflam = getLamNfTwoLoop(q);
    val = alphaS(q, nflam.second, nflam.first);
  }
  // special handling if the scale is less than Qmin
  else {
    switch (_asType) {
    case 1: 
      // flat, zero; the default type with no NP effects.
      val = 0.; 
      break;
    case 2: 
      // flat, non-zero alpha_s = alpha_s(q2min).
      val = _val0;
      break;
    case 3: 
      // linear in q
      val = _val0*q/_qmin;
      break;
    case 4:
      // quadratic in q
      val = _val0*sqr(q/_qmin);
      break;
    case 5:
      // quadratic in q, starting off at asMaxNP, ending on as(qmin)
      val = (_val0 - _asMaxNP)*sqr(q/_qmin) + _asMaxNP;
      break;
    case 6:
      // just asMaxNP and constant
      val = _asMaxNP;
      break;
    }
  }
  return val;
}

double ShowerAlphaQCD::overestimateValue() const {
  return _alphamin;
}

double ShowerAlphaQCD::ratio(const Energy2 scale, double factor) const {
  Energy q = scaleFactor()*factor*sqrt(scale);
  double val(0.);
  // normal case
  if (q >= _qmin) {
    pair<short,Energy> nflam = getLamNfTwoLoop(q);
    val = alphaS(q, nflam.second, nflam.first);
  }
  // special handling if the scale is less than Qmin
  else {
    switch (_asType) {
    case 1:
      // flat, zero; the default type with no NP effects.
      val = 0.; 
      break;
    case 2:
      // flat, non-zero alpha_s = alpha_s(q2min).
      val = _val0;
      break;
    case 3:
      // linear in q
      val = _val0*q/_qmin;
      break;
    case 4:
      // quadratic in q
      val = _val0*sqr(q/_qmin);
      break;
    case 5:
      // quadratic in q, starting off at asMaxNP, ending on as(qmin)
      val = (_val0 - _asMaxNP)*sqr(q/_qmin) + _asMaxNP;
      break;
    case 6:
      // just asMaxNP and constant
      val = _asMaxNP;
      break;
    }
  }
  // denominator 
  return val/_alphamin;
}

string ShowerAlphaQCD::value (string scale) {
  istringstream readscale(scale);
  double inScale; readscale >> inScale;
  Energy theScale = inScale * GeV;
  initialize();
  ostringstream showvalue ("");
  showvalue << "alpha_s (" << theScale/GeV << " GeV) = "
	    << value (sqr(theScale));
  return showvalue.str();
}

pair<short, Energy> ShowerAlphaQCD::getLamNfTwoLoop(Energy q) const {
  short nf = 6;
  // get lambda and nf according to the thresholds
  if      (q < _thresholds[1]) nf = 3;
  else if (q < _thresholds[2]) nf = 4;
  else if (q < _thresholds[3]) nf = 5;
  return pair<short,Energy>(nf, _lambda[nf-3]);
}

Energy ShowerAlphaQCD::computeLambda(Energy match,
				     double alpha,
				     unsigned int nflav) const {
  Energy lamtest=200.0*MeV;
  double xtest;
  unsigned int ntry=0;
  do {
    ++ntry;
    xtest  = log(sqr(match/lamtest));
    xtest += (alpha-alphaS(match,lamtest,nflav))/derivativealphaS(match,lamtest,nflav);
    Energy newLambda = match/exp(0.5*xtest);
    lamtest = newLambda<match ? newLambda : 0.5*(lamtest+match);
  }
  while(abs(alpha-alphaS(match,lamtest,nflav)) > _tolerance && ntry < _maxtry);
  return lamtest;
}

double ShowerAlphaQCD::derivativealphaS(Energy q, Energy lam, int nf) const {
  using Constants::pi;
  double lx = log(sqr(q/lam));
  double b0 = 11. - 2./3.*nf;
  double b1 = 51. - 19./3.*nf;
  double b2 = 2857. - 5033./9.*nf + 325./27.*sqr(nf);
  if(_nloop==1)
    return -4.*pi/(b0*sqr(lx));
  else if(_nloop==2)
    return -4.*pi/(b0*sqr(lx))*(1.+2.*b1/sqr(b0)/lx*(1.-2.*log(lx)));
  else
    return -4.*pi/(b0*sqr(lx))*
      (1.  + 2.*b1/sqr(b0)/lx*(1.-2.*log(lx))
       + 4.*sqr(b1)/(sqr(sqr(b0))*sqr(lx))*(1. - 2.*log(lx)
					    + 3.*(sqr(log(lx) - 0.5)+b2*b0/(8.*sqr(b1))-1.25)));
}

double ShowerAlphaQCD::alphaS(Energy q, Energy lam, int nf) const {
  using Constants::pi;
  double lx(log(sqr(q/lam)));
  double b0 = 11. - 2./3.*nf;
  double b1 = 51. - 19./3.*nf;
  double b2 = 2857. - 5033./9.*nf + 325./27.*sqr(nf);
  // one loop
  if(_nloop==1)
    {return 4.*pi/(b0*lx);}
  // two loop
  else if(_nloop==2) {
    return 4.*pi/(b0*lx)*(1.-2.*b1/sqr(b0)*log(lx)/lx);
  }
  // three loop
  else
    {return 4.*pi/(b0*lx)*(1.-2.*b1/sqr(b0)*log(lx)/lx + 
			   4.*sqr(b1)/(sqr(sqr(b0))*sqr(lx))*
			   (sqr(log(lx) - 0.5) + b2*b0/(8.*sqr(b1)) - 5./4.));}
}

