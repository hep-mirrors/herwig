// -*- C++ -*-
//
// ShowerAlphaQCDNP.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQCDNP class.
//
#include "ShowerAlphaQCDNP.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Config/Constants.h"
#include "Herwig/Utilities/AlphaS.h"
#include "gsl/gsl_sf_lambert.h"

using namespace Herwig;
using Herwig::Math::alphaS;
using Herwig::Math::derivativeAlphaS;

DescribeClass<ShowerAlphaQCDNP,ShowerAlpha>
describeShowerAlphaQCDNP("Herwig::ShowerAlphaQCDNP","HwShower.so");

IBPtr ShowerAlphaQCDNP::clone() const {
  return new_ptr(*this);
}

IBPtr ShowerAlphaQCDNP::fullclone() const {
  return new_ptr(*this);
}

void ShowerAlphaQCDNP::persistentOutput(PersistentOStream & os) const {
  os << ounit(_qmin,GeV) << _asMaxNP << _nfMaxNP
     << ounit(_thresholds,GeV) << ounit(_lambda,GeV)
     << _nloop << _thresopt 
     << _alphain << _tolerance << _maxtry
     << _npPower << ounit(_optInputScale,GeV) << _quarkBranching
     << ounit(_quarkMasses,GeV);
}

void ShowerAlphaQCDNP::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_qmin,GeV) >> _asMaxNP >> _nfMaxNP
     >> iunit(_thresholds,GeV) >> iunit(_lambda,GeV)
     >> _nloop >> _thresopt 
     >> _alphain >> _tolerance >> _maxtry
     >> _npPower >> iunit(_optInputScale,GeV) >> _quarkBranching
     >> iunit(_quarkMasses,GeV);
}

void ShowerAlphaQCDNP::Init() {

  static ClassDocumentation<ShowerAlphaQCDNP> documentation
    ("This (concrete) class describes the QCD alpha running.");

  // default such that as(qmin) = 1 in the current parametrization.
  static Parameter<ShowerAlphaQCDNP,Energy> intQmin
    ("Qmin", "Q < Qmin is treated with NP parametrization as of (unit [GeV])",
     &ShowerAlphaQCDNP::_qmin, GeV, 1*GeV, 1*GeV,
     100.0*GeV,false,false,false);

  static Parameter<ShowerAlphaQCDNP,unsigned int> interfaceNumberOfLoops
    ("NumberOfLoops",
     "The number of loops to use in the alpha_S calculation",
     &ShowerAlphaQCDNP::_nloop, 3, 1, 3,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCDNP,double> interfaceAlphaMZ
    ("AlphaIn",
     "The input value of the strong coupling at the chosen InputScale (default: MZ)",
     &ShowerAlphaQCDNP::_alphain, 0.118, 0.1, 0.2,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCDNP,double> interfaceTolerance
    ("Tolerance",
     "The tolerance for discontinuities in alphaS at thresholds.",
     &ShowerAlphaQCDNP::_tolerance, 1e-10, 1e-20, 1e-4,
     false, false, Interface::limited);

  static Parameter<ShowerAlphaQCDNP,unsigned int> interfaceMaximumIterations
    ("MaximumIterations",
     "The maximum number of iterations for the Newton-Raphson method to converge.",
     &ShowerAlphaQCDNP::_maxtry, 100, 10, 1000,
     false, false, Interface::limited);

  static Switch<ShowerAlphaQCDNP,bool> interfaceThresholdOption
    ("ThresholdOption",
     "Whether to use the consistuent or normal masses for the thresholds",
     &ShowerAlphaQCDNP::_thresopt, true, false, false);
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

  static Command<ShowerAlphaQCDNP> interfaceValue
    ("Value",
     "",
     &ShowerAlphaQCDNP::value, false);

  static Command<ShowerAlphaQCDNP> interfacecheck
    ("check",
     "check",
     &ShowerAlphaQCDNP::check, false);

  static Command<ShowerAlphaQCDNP> interfacecheckNf
    ("checkNf",
     "checkNf",
     &ShowerAlphaQCDNP::checkNf, false);  

    static Command<ShowerAlphaQCDNP> interfacecheckscaleNP
    ("checkscaleNP",
     "checkscaleNP",
     &ShowerAlphaQCDNP::checkscaleNP, false);  

  static Parameter<ShowerAlphaQCDNP,Energy> interfaceInputScale
    ("InputScale",
     "An optional input scale. MZ will be used if not set.",
     &ShowerAlphaQCDNP::_optInputScale, GeV, 91.1876_GeV, ZERO, 0*GeV,
     false, false, Interface::lowerlim);

  static ParVector<ShowerAlphaQCDNP,Energy> interfaceQuarkMasses
    ("QuarkMasses",
     "The quark masses to be used instead of the masses set in the particle data.",
     &ShowerAlphaQCDNP::_quarkMasses, GeV, -1, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  
  static Switch<ShowerAlphaQCDNP,bool> interfaceQuarkBranching
    ("QuarkBranching",
     "True, if this coupling is used in a gluon to qqbar branching.",
     &ShowerAlphaQCDNP::_quarkBranching, false, false, false);
  static SwitchOption interfaceQuarkBranchingYes
    (interfaceQuarkBranching,
     "Yes",
     "Use in gluon to qqbar branching.",
     true);
  static SwitchOption interfaceQuarkBranchingNo
    (interfaceQuarkBranching,
     "No",
     "Use in gluon emission.",
     false);
  
  static Parameter<ShowerAlphaQCDNP,double> interfaceNPPower
    ("NPPower",
     "The non-perturbative power law",
     &ShowerAlphaQCDNP::_npPower, 2.0, 1.0, 0,
     false, false, Interface::lowerlim);

}

void ShowerAlphaQCDNP::doinit() {
  ShowerAlpha::doinit();
  // calculate the value of 5-flavour lambda 
  // evaluate the initial
  // value of Lambda from alphas if needed using Newton-Raphson
  _lambda[2]=computeLambda(_optInputScale,_alphain,5);  

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
  _lambda[3]=computeLambda(_thresholds[3],alphaS(_thresholds[3],_lambda[2],5,_nloop),6);
  // bottom threshold
  // compute 4 flavour lambda by matching at bottom mass using Newton Raphson
  _lambda[1]=computeLambda(_thresholds[2],alphaS(_thresholds[2],_lambda[2],5,_nloop),4);
  // charm threshold
  // compute 3 flavour lambda by matching at charm mass using Newton Raphson
  _lambda[0]=computeLambda(_thresholds[1],alphaS(_thresholds[1],_lambda[1],4,_nloop),3);

  //Energy q = scaleNPMin();
  //auto nflam = getLamNf(q);
  //_asMaxNP = alphaS(q, nflam.second, nflam.first, _nloop);

  Energy q = scaleThatmMinimizesScaleNP() ;
  _asMaxNP = value(q*q)*1.05;

  //_absoluteCutoff=_qmin*pow(100.*log(10.),-1./_npPower);   
  _absoluteCutoff=2.0*GeV;

  _nfMaxNP = nfNP(sqr(_absoluteCutoff));

}

string ShowerAlphaQCDNP::check(string args) {

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

    //out << (Q/GeV) << " " << value(Q*Q) << " " << _asMaxNP << "\n";
    out << (Q/GeV) << " " << value(Q*Q)  << "\n";


  }

  return "alpha_s check finished";

}

string ShowerAlphaQCDNP::checkNf(string args) {

  doinit();

  istringstream argin(args);

  double Q_low, Q_high;
  long n_steps;

  argin >> Q_low >> Q_high >> n_steps;

  string fname;
  argin >> fname;

  Repository::clog() << "checking nfNP in range [" << Q_low << "," << Q_high << "] GeV in "
		     << n_steps << " steps.\nResults are written to " << fname << "\n";

  double step_width = (Q_high-Q_low)/n_steps;

  ofstream out (fname.c_str());

  for (long k = 0; k <= n_steps; ++k) {

    Energy Q = Q_low*GeV + k*step_width*GeV;

    //out << (Q/GeV) << " " << nfNP(Q*Q) << " " << _nfMaxNP << "\n";
    out << (Q/GeV) << " " << nfNP(Q*Q) << "\n";

  }

  return "nfNP check finished";

}



string ShowerAlphaQCDNP::checkscaleNP(string args) {

  doinit();

  istringstream argin(args);

  double Q_low, Q_high;
  long n_steps;

  argin >> Q_low >> Q_high >> n_steps;

  string fname;
  argin >> fname;

  Repository::clog() << "checking scaleNP (in GeV) in range [" << Q_low << "," << Q_high << "] GeV in "
         << n_steps << " steps.\nResults are written to " << fname << "\n";

  double step_width = (Q_high-Q_low)/n_steps;

  ofstream out (fname.c_str());

  for (long k = 0; k <= n_steps; ++k) {

    Energy Q = Q_low*GeV + k*step_width*GeV;

    out << (Q/GeV) << " " << scaleNP(scaleFactor()*Q)/(1.*GeV) << "\n";

  }

  return "scaleNP check finished";

}


double ShowerAlphaQCDNP::value(const Energy2 scale) const {

  Energy q = scaleFactor()*sqrt(scale);
  
  if ( q > _absoluteCutoff ) {  
    auto nflam = getLamNf(q);
    return alphaS(scaleNP(q), nflam.second, nflam.first, _nloop);    
  }


  return 0.;

}


string ShowerAlphaQCDNP::value(string scale) {
  istringstream readscale(scale);
  double inScale; readscale >> inScale;
  Energy theScale = inScale * GeV;
  initialize();
  ostringstream showvalue ("");
  showvalue << "alpha_s (" << theScale/GeV << " GeV) = "
	    << value (sqr(theScale));
  return showvalue.str();
}

pair<short, Energy> ShowerAlphaQCDNP::getLamNf(Energy q) const {
  short nf = 6;
  // get lambda and nf according to the thresholds
  if      (q < _thresholds[1]) nf = 3;
  else if (q < _thresholds[2]) nf = 4;
  else if (q < _thresholds[3]) nf = 5;
  return pair<short,Energy>(nf, _lambda[nf-3]);
}

Energy ShowerAlphaQCDNP::computeLambda(Energy match,
				     double alpha,
				     unsigned int nflav) const {
	
  Energy lamtest=200.0*MeV;
  double xtest;
  unsigned int ntry=0;
  do {
    ++ntry;
    xtest  = log(sqr(match/lamtest));
    xtest += (alpha-alphaS(match,lamtest,nflav,_nloop))/derivativeAlphaS(match,lamtest,nflav,_nloop);
    Energy newLambda = match/exp(0.5 *xtest);
    lamtest = newLambda<match ? newLambda : 0.5*(lamtest+match);
  }
  while(abs(alpha-alphaS(match,lamtest,nflav,_nloop)) > _tolerance && ntry < _maxtry);
  return lamtest;

}

// --------------------------------------------------------------------------------
// main modification for nonperturbative running in terms of changed argument and
// nonperturbative number of flavours
// --------------------------------------------------------------------------------

double ShowerAlphaQCDNP::nfNP(const Energy2 q2) const {
  Energy q=sqrt(q2);

  if (q<_absoluteCutoff){
    return 0.;
  }
  else{

  //get p-nf
  double nf=getLamNf(q).first;

  //beta function coefficients; b=sum_i (alphaS/(4Pi))^i * sum_j b_ij *nf^j
  using Constants::pi;
  double b00 = 11.;
  double b01 = - 2./3.;
  double b10 = 51.;
  double b11 = - 19./3.;
  double c = (b00+b10*value(q2)/(4.*pi))/(b01+b11*value(q2)/(4.*pi));

  //absolute number of np-flavors
  double nfNPabs=(c+nf)*(q*scaleNPDerivative(q)/scaleNP(q))-c;

  //return ratio of np-flavors over p-flavors
  return nfNPabs/nf;
  }
}

Energy ShowerAlphaQCDNP::scaleNP(const Energy q) const {
  return q*(1.+(_qmin/q)*(exp(pow(_qmin/q,_npPower))-1.));
}

double ShowerAlphaQCDNP::scaleNPDerivative(const Energy q) const {
  return 1.-exp(pow(_qmin/q,_npPower))*_npPower*pow(_qmin/q,1.+_npPower);
}

Energy ShowerAlphaQCDNP::scaleThatmMinimizesScaleNP() const {
  double x=(1./_npPower)*gsl_sf_lambert_W0(pow(_npPower,1./(1.+_npPower))/(1.+_npPower));
  return _qmin*exp(x)*pow(_npPower,1./(1.+_npPower));

}

Energy ShowerAlphaQCDNP::scaleNPMin() const {
  double x=((1.+_npPower)/_npPower)*gsl_sf_lambert_W0(pow(_npPower,1./(1.+_npPower))/(1.+_npPower));
  return _qmin*(-1.+exp(x)+pow(x,-1./_npPower));
}
