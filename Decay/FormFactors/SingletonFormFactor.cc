// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SingletonFormFactor class.
//

#include "SingletonFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void SingletonFormFactor::doinit() throw(InitException) {
  BaryonFormFactor::doinit();
  if(numberOfFactors()!=_polemass.size())
    throw InitException() << "Inconsistent parameters in SingletonFormFactor::doinit()"
			  << Exception::abortnow;
  // calculate the constants for the form-factors
  int id0,id1;
  _xi.resize(0);_NmM.resize(0);_mquark.resize(0);
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    // id codes for the particles
    particleID(ix,id0,id1);
    if((abs(id0)==5122&&abs(id1)==4122)||(abs(id0)==5132&&abs(id1)==4132)||
       (abs(id0)==5232&&abs(id1)==4232)) {
      _mquark.push_back(_mcharm);
      _xi.push_back(1.);
      _NmM.push_back(1.);
    }
    else if((abs(id0)==5222&&abs(id1)==4222)||(abs(id0)==5212&&abs(id1)==4212)||
	    (abs(id0)==5112&&abs(id1)==4112)||(abs(id0)==5332&&abs(id1)==4332)||
	    (abs(id0)==5312&&abs(id1)==4312)||(abs(id0)==5322&&abs(id1)==4322)) {
      _mquark.push_back(_mcharm);
      _xi.push_back(-1./3.);
      _NmM.push_back(1.);
    }
    else if(abs(id0)==4122&&abs(id1)==3122) {
      _mquark.push_back(_mstrange);
      _xi.push_back(1.);
      _NmM.push_back(sqrt(2./3.)*sin(_thetalambda));
    }
    else if((abs(id0)==4222&&abs(id1)==3222)||(abs(id0)==4212&&abs(id1)==3212)||
	    (abs(id0)==4112&&abs(id1)==3112)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(-1./3.);
      _NmM.push_back(sqrt(2./3.)*cos(_thetasigma));
    }
    else if((abs(id0)==4132&&abs(id1)==3312)||(abs(id0)==4232&&abs(id1)==3322)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(1.);
      _NmM.push_back(1./sqrt(2.)*sin(_thetaxi));
    }
    else if((abs(id0)==4312&&abs(id1)==3322)||(abs(id0)==4322&&abs(id1)==3312)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(-1./3.);
      _NmM.push_back(1./sqrt(6.)*cos(_thetaxip));
    }
    else
      throw InitException() << "Unknown mode in SingletonFormFactor::doinit()"
			    << Exception::abortnow;
  }
}

void SingletonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _mcharm << _mstrange <<  _thetalambda << _thetasigma << _thetaxi 
     << _thetaxip << _polemass << _xi << _NmM << _mquark;
}

void SingletonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _mcharm >> _mstrange >>  _thetalambda >> _thetasigma >> _thetaxi 
     >> _thetaxip >> _polemass >> _xi >> _NmM >> _mquark;
}

ClassDescription<SingletonFormFactor> SingletonFormFactor::initSingletonFormFactor;
// Definition of the static class description member.

void SingletonFormFactor::Init() {

  static ClassDocumentation<SingletonFormFactor> documentation
    ("The SingletonFormFactor class implements the"
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
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int, int, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a)
{
  InvEnergy ratio(0.5/m0);
  // all factors divided by sqrt(4.*m0*m1) normalisation  
  double gbar(m1*_xi[iloc]*_NmM[iloc]/_mquark[iloc]);
  double abar(_xi[iloc]*_NmM[iloc]);
  InvEnergy apbar(-ratio*(m1/_mquark[iloc]-1.)*_xi[iloc]*_NmM[iloc]);
  InvEnergy ambar(apbar);
  InvEnergy gpbar(-ratio*_NmM[iloc]*((_xi[iloc]*m1/_mquark[iloc]-1.)
				  +0.5*(m1-m0)/_mquark[iloc]*(1.-_xi[iloc])));
  InvEnergy gmbar(-ratio*_NmM[iloc]*((_xi[iloc]*m1/_mquark[iloc]-1.)
				  +0.5*(m1+m0)/_mquark[iloc]*(1.-_xi[iloc])));
  // energy dependence
  double ymax((1.-m1/m0));        ymax *=ymax;
  double yres(_polemass[iloc]/m0);yres *=yres;
  double y(q2/m0/m0),efact((ymax-yres)/(y-yres));
  f1v = efact*(gbar+(m0+m1)*gpbar);
  f2v = efact*gpbar*(m0+m1);
  f3v = efact*gmbar*(m0+m1);
  f1a = efact*(-abar+(m0-m1)*apbar);
  f2a =-efact*apbar*(m0+m1);
  f3a =-efact*ambar*(m0+m1);
}

  void SingletonFormFactor::dataBaseOutput(ofstream & output,bool header,
					   bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create /Herwig++/SingletonFormFactor " << fullName() << " \n ";}
  output << "set " << fullName() << ":CharmMass " << _mcharm/GeV << " \n";
  output << "set " << fullName() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "set " << fullName() << ":ThetaLambda " << _thetalambda << " \n";
  output << "set " << fullName() << ":ThetaSigma " << _thetasigma << " \n";
  output << "set " << fullName() << ":ThetaXi " << _thetaxi << " \n";
  output << "set " << fullName() << ":ThetaXiPrime " << _thetaxip << " \n";
  for(unsigned int ix=0;ix<_polemass.size();++ix)
    {
      if(ix<initialModes())
	{output << "set " << fullName() << ":PoleMass " << ix << "  " 
		<< _polemass[ix]/GeV << endl;}
      else
	{output << "insert " << fullName() << ":PoleMass "<< ix << "  " 
		<< _polemass[ix]/GeV << endl;}
    }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
