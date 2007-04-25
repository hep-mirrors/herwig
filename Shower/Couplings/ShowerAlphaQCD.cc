// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQCD class.
//
#include "ShowerAlphaQCD.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

void ShowerAlphaQCD::persistentOutput(PersistentOStream & os) const {
  os << _asType << _asMaxNP << _qmin << _nloop << _lambdaopt << _thresopt 
     << _lambdain << _alphain << _inopt
     << _tolerance << _maxtry << _alphamin;
  for(unsigned int ix=0;ix<4;++ix)
    os << _thresholds[ix] << _lambda[ix];
}

void ShowerAlphaQCD::persistentInput(PersistentIStream & is, int) {
  is >> _asType >> _asMaxNP >> _qmin >> _nloop >> _lambdaopt >> _thresopt
     >> _lambdain >> _alphain >> _inopt
     >> _tolerance >> _maxtry >> _alphamin;
  for(unsigned int ix=0;ix<4;++ix)
    is >> _thresholds[ix] >> _lambda[ix];
}

ClassDescription<ShowerAlphaQCD> ShowerAlphaQCD::initShowerAlphaQCD;
// Definition of the static class description member.

void ShowerAlphaQCD::Init() {

  static ClassDocumentation<ShowerAlphaQCD> documentation
    ("This (concrete) class describes the QCD alpha running.");

  static Switch<ShowerAlphaQCD, int> intAsType
    ("NPAlphaS",
     "Behaviour of AlphaS in the NP region",
     &ShowerAlphaQCD::_asType, 1, false, false);
  static SwitchOption intAsTypeZero
    (intAsType, "AsTypeZero","zero below Q_min", 1);
  static SwitchOption intAsTypeConst
    (intAsType, "AsTypeConst","const as(qmin) below Q_min", 2);
  static SwitchOption intAsTypeLin
    (intAsType, "AsTypeLin ","growing linearly below Q_min", 3);
  static SwitchOption intAsTypeQuad
    (intAsType, "AsTypeQuad","growing quadratically below Q_min", 4);
  static SwitchOption intAsTypeExx1
    (intAsType, "AsTypeExx1 ", "quad from AlphaMaxNP down to as(Q_min)", 5);
  static SwitchOption intAsTypeExx2
    (intAsType, "AsTypeExx2 ", "const = AlphaMaxNP below Q_min", 6);

  // default such that as(qmin) = 1 in the current parametrization.
  // min = Lambda3
  static Parameter<ShowerAlphaQCD,Energy> intQmin
    ("Qmin", "Q < Qmin is treated with NP parametrization as of (unit [GeV])",
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

}

void ShowerAlphaQCD::doinit() throw(InitException) {
  ShowerAlpha::doinit();
  // calculate the value of 5-flavour lambda 
  // evaluate the initial
  // value of Lambda from alphas if needed using Newton-Raphson
  if(_inopt)
    {_lambda[2]=computeLambda(getParticleData(ParticleID::Z0)->mass(),_alphain,5);}
  // otherwise it was an input parameter
  else{_lambda[2]=_lambdain;}
  // convert lambda to the Monte Carlo scheme if needed
  if(_lambdaopt){_lambda[2] *=exp(0.5*(67.-3.*sqr(pi)-50./3.)/23.)/sqrt(2.);}
  // compute the threshold matching
  // top threshold
  for(int ix=1;ix<4;++ix) {
    if(_thresopt)
      _thresholds[ix]=getParticleData(ix+3)->mass();
    else
      _thresholds[ix]=getParticleData(ix+3)->constituentMass();
  }
  // compute 6 flavour lambda by matching at top mass using Newton Raphson
  _lambda[3]=computeLambda(_thresholds[3],alphaS(_thresholds[3],_lambda[2],5),6);
  // bottom threshold
  // compute 4 flavour lambda by matching at bottom mass using Newton Raphson
  _lambda[1]=computeLambda(_thresholds[2],alphaS(_thresholds[2],_lambda[2],5),4);
  // charm threshold
  // compute 3 flavour lambda by matching at charm mass using Newton Raphson
  _lambda[0]=computeLambda(_thresholds[1],alphaS(_thresholds[1],_lambda[1],4),3);
  // final threshold is qmin
  _thresholds[0]=_qmin;
  // compute the maximum value of as 
  if ( _asType < 5 ) _alphamin = value(sqr(_qmin)); // gives as = 1
  else _alphamin = _asMaxNP; 
  // check consistency lambda_3 < qmin
  if(_lambda[0]>_qmin)
    Throw<InitException>() << "The value of Qmin is less than Lambda_3 in"
			   << " ShowerAlphaQCD::doinit " << Exception::abortnow;
}

double ShowerAlphaQCD::value(const Energy2 scale) const {
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
  cout << "q = " << q/GeV << "GeV, qmin = " << _qmin/GeV << " GeV, value = " << scaleFactor()*val << endl; 
  return scaleFactor() * val;
}

double ShowerAlphaQCD::overestimateValue() const {
  cout << "overestimateValue = " << scaleFactor()*_alphamin << endl; 
  return scaleFactor() * _alphamin; 
}

double ShowerAlphaQCD::ratio(const Energy2 scale) const {
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
  cout << "q = " << q/GeV << "GeV, qmin = " << _qmin/GeV << "GeV, ratio = " << val << "/" << _alphamin << " = " << val/_alphamin << endl; 
  return val/_alphamin;  
}
