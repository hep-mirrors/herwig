// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SingletonFormFactor class.
//

#include "SingletonFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SingletonFormFactor.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;

SingletonFormFactor::~SingletonFormFactor() {}

void SingletonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _mcharm << _mstrange <<  _thetalambda << _thetasigma << _thetaxi 
     << _thetaxip << _polemass << _abar << _gbar << _apbar << _ambar 
     << _gpbar << _gmbar;
}

void SingletonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _mcharm >> _mstrange >>  _thetalambda >> _thetasigma >> _thetaxi 
     >> _thetaxip >> _polemass >> _abar >> _gbar >> _apbar >> _ambar 
     >> _gpbar >> _gmbar;
}

ClassDescription<SingletonFormFactor> SingletonFormFactor::initSingletonFormFactor;
// Definition of the static class description member.

void SingletonFormFactor::Init() {

  static ClassDocumentation<SingletonFormFactor> documentation
    ("The \\classname{SingletonFormFactor} class implements the"
     " form-factors of PRD43, 2939 for the decay of spin-1/2 baryons"
     " containing one heavy quark.");

  static Parameter<SingletonFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark",
     &SingletonFormFactor::_mcharm, GeV, 1.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<SingletonFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark",
     &SingletonFormFactor::_mstrange, GeV, 0.51*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaLambda
    ("ThetaLambda",
     "The mixing angle for the Lambda",
     &SingletonFormFactor::_thetalambda, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaSigma
    ("ThetaSigma",
     "The mixing angle for the Sigma",
     &SingletonFormFactor::_thetasigma, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaXi
    ("ThetaXi",
     "The mixing angle for the Xi",
     &SingletonFormFactor::_thetaxi, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaXiPrime
    ("ThetaXiPrime",
     "The mixing angle for the Xi'",
     &SingletonFormFactor::_thetaxip, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static ParVector<SingletonFormFactor,Energy> interfacePoleMass
    ("PoleMass",
     "The mass for the energy dependence of the form-factors.",
     &SingletonFormFactor::_polemass,
     1.*GeV, 0, 0, -10.*GeV, 10.*GeV, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void SingletonFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int id0, int id1, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  // energy dependence
  double ymax = (1.-m1/m0); ymax *=ymax;
  double yres = _polemass[iloc]/m0;yres*=yres;
  double y    = q2/m0/m0;
  InvEnergy efact=(ymax-yres)/(y-yres)/sqrt(4.*m0*m1);
  f1v = efact*(_gbar[iloc]+(m0+m1)*_gpbar[iloc]);
  f2v = efact*_gpbar[iloc]*(m0+m1);
  f3v = efact*_gmbar[iloc]*(m0+m1);
  f1a = efact*(-_abar[iloc]+(m0-m1)*_apbar[iloc]);
  f2a =-efact*_apbar[iloc]*(m0+m1);
  f3a =-efact*_ambar[iloc]*(m0+m1);
}

void SingletonFormFactor::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/SingletonFormFactor " << fullName() << " \n ";
  output << "set " << fullName() << ":CharmMass " << _mcharm/GeV << " \n";
  output << "set " << fullName() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "set " << fullName() << ":ThetaLambda " << _thetalambda << " \n";
  output << "set " << fullName() << ":ThetaSigma " << _thetasigma << " \n";
  output << "set " << fullName() << ":ThetaXi " << _thetaxi << " \n";
  output << "set " << fullName() << ":ThetaXiPrime " << _thetaxip << " \n";
  for(unsigned int ix=0;ix<_polemass.size();++ix)
    {
      if(ix<13){output << "set " << fullName() << ":PoleMass " << ix << "  " 
		       << _polemass[ix]/GeV << endl;}
      else{output << "insert " << fullName() << ":PoleMass "<< ix << "  " 
		  << _polemass[ix]/GeV << endl;}
    }
}
  
}
