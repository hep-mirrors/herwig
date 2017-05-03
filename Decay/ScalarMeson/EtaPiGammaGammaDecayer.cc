// -*- C++ -*-
//
// EtaPiGammaGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiGammaGammaDecayer class.
//
#include "EtaPiGammaGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void EtaPiGammaGammaDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _etamax  = mode(0)->maxWeight();
    _etapmax = mode(1)->maxWeight();
  }
}

EtaPiGammaGammaDecayer::EtaPiGammaGammaDecayer()
  : _grhoomega(12.924/GeV), _fpi(130.7*MeV),_rhomass(771.1*MeV),
    _rhowidth(149.2*MeV),_grho(_rhomass/_fpi),_mpi(ZERO),_rhoconst(0.),
    _localparameters(true),_ratiofpif8(1./1.3),_ratiofpif0(1./1.04),
    _theta(-Constants::pi/9.),_etamax(2.36858),_etapmax(0.006),
    _dconst(2), _econst(2) {
  // intermediates
  generateIntermediates(false);
}

void EtaPiGammaGammaDecayer::doinit() {
  DecayIntegrator::doinit();
  // set rho parameters if needed
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters) {
    _rhomass  = rho->mass();
    _rhowidth = rho->width();
  }
  // constant for the running rho width
  _mpi=getParticleData(ParticleID::pi0)->mass();
  Energy pcm =Kinematics::pstarTwoBodyDecay(_rhomass,_mpi,_mpi);
  _rhoconst=_rhomass*_rhomass*_rhowidth/(pcm*pcm*pcm);
  // set the prefactors
  double conv(sqrt(4.*Constants::pi*SM().alphaEM()));
  conv *=_fpi*_fpi*_grho/_rhomass/_rhomass;
  InvEnergy2 pre(2.*sqrt(3.)/9.*sqr(_grhoomega*conv));
  double fact[2];
  // constants for eta
  fact[0] = _ratiofpif8*cos(_theta)-sqrt(2.)*_ratiofpif0*sin(_theta);
  // constants for eta'
  fact[1] = _ratiofpif8*sin(_theta)+sqrt(2.)*_ratiofpif0*cos(_theta);
  for(unsigned int ix=0;ix<2;++ix) {
    _dconst[ix]=fact[ix]*pre;
    _econst[ix]=fact[ix]*pre;
  }
  // set up the phsae space for the decays
  tPDPtr eta[2]={getParticleData(ParticleID::eta),getParticleData(ParticleID::etaprime)};
  tPDVector extpart;extpart.resize(4);
  extpart[1] = getParticleData(ParticleID::pi0);
  extpart[2] = getParticleData(ParticleID::gamma);
  extpart[3] = getParticleData(ParticleID::gamma);
  vector<double> dummyweights(2,0.5);
  DecayPhaseSpaceChannelPtr newchannel;
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<2;++ix) {
    extpart[0] = eta[ix];
    mode = new_ptr(DecayPhaseSpaceMode(extpart,this));
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,2);
    newchannel->addIntermediate(rho,0,0.0, 1,3);
    mode->addChannel(newchannel);
    newchannel=new_ptr(DecayPhaseSpaceChannel(mode));
    newchannel->addIntermediate(extpart[0],0, 0.0,-1,3);
    newchannel->addIntermediate(rho,0,0.0, 1,2);
    mode->addChannel(newchannel);
    addMode(mode, (ix==0 ? _etamax : _etapmax),dummyweights);
  }
}

int EtaPiGammaGammaDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  cc=false;
  int id;
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npi0(0),ngamma(0);
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0)         ++npi0;
    else if(id==ParticleID::gamma)  ++ngamma;
  }
  if(!(npi0==1&&ngamma==2)) return -1;
  // number of the mode
  switch (parent->id()) {
  case ParticleID::eta     : return 0;
  case ParticleID::etaprime: return 1;
  default: return -1;
  }
}

void EtaPiGammaGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_grhoomega,1/GeV)<< ounit(_fpi,GeV)<< _grho 
     << ounit(_rhomass,GeV)<< ounit(_rhowidth,GeV)<< _localparameters 
     << _ratiofpif8 << _ratiofpif0 << _theta << _etamax << _etapmax 
     << _rhoconst << ounit(_mpi,GeV) << ounit(_dconst,1/GeV2) 
     << ounit(_econst,1/GeV2);
}

void EtaPiGammaGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_grhoomega,1/GeV) >> iunit(_fpi,GeV)>> _grho 
     >> iunit(_rhomass,GeV)>> iunit(_rhowidth,GeV)>> _localparameters 
     >> _ratiofpif8 >> _ratiofpif0 >> _theta >> _etamax >> _etapmax 
     >> _rhoconst >> iunit(_mpi,GeV) >> iunit(_dconst,1/GeV2) 
     >> iunit(_econst,1/GeV2);
}

ClassDescription<EtaPiGammaGammaDecayer> 
EtaPiGammaGammaDecayer::initEtaPiGammaGammaDecayer;
// Definition of the static class description member.

void EtaPiGammaGammaDecayer::Init() {

  static ClassDocumentation<EtaPiGammaGammaDecayer> documentation
    ("The EtaPiGammaGammaDecayer class implements a VMD model for the"
     " decay of the eta or etaprime to a pion and two photons.",
     "The decays of $\\eta,\\eta'\\to\\pi^0\\gamma\\gamma$ were simulated using"
     " the matrix elements of \\cite{Holstein:2001bt}",
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Parameter<EtaPiGammaGammaDecayer,InvEnergy> interfacegrhoomega
    ("grhoomega",
     "The couping of the rho, omega and a pion",
     &EtaPiGammaGammaDecayer::_grhoomega, 1./GeV, 12.924/GeV, ZERO, 100./GeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &EtaPiGammaGammaDecayer::_fpi, MeV, 130.7*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfacegrho
    ("grho",
     "Rho decay constant",
     &EtaPiGammaGammaDecayer::_grho, 5.9, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho meson",
     &EtaPiGammaGammaDecayer::_rhomass, MeV, 771.1*MeV, 500.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho meson",
     &EtaPiGammaGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF8
    ("RatioFpiF8",
     "The ratio of the decay constant Fpi to F8",
     &EtaPiGammaGammaDecayer::_ratiofpif8, 1./1.3, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF0
    ("RatioFpiF0",
     "The ratio of the decay constant Fpi to F0",
     &EtaPiGammaGammaDecayer::_ratiofpif0, 1./1.04, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceTheta
    ("Theta",
     "The eta etaprime mixing angle",
     &EtaPiGammaGammaDecayer::_theta, -Constants::pi/9., -Constants::pi, Constants::pi,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaMax
    ("EtaMax",
     "THe maximum weight for the eta decay",
     &EtaPiGammaGammaDecayer::_etamax, 1.35, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaPrimeMax
    ("EtaPrimeMax",
     "THe maximum weight for the eta prime decay",
     &EtaPiGammaGammaDecayer::_etapmax, 0.006, -1.0e12, 1.0e12,
     false, false, false);

  static Switch<EtaPiGammaGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the parameters",
     &EtaPiGammaGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

double EtaPiGammaGammaDecayer::me2(const int,const Particle & inpart,
				   const ParticleVector& decay,
				   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix+1],
					    outgoing,true,true);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::calculateWaveFunctions(_vectors[ix],decay[ix+1],
					       outgoing,true);
  // dot products we need
  Energy2 q1dotq2(decay[1]->momentum()*decay[2]->momentum()),
    pdotq1(inpart.momentum()*decay[1]->momentum()),
    pdotq2(inpart.momentum()*decay[2]->momentum());
  complex<Energy> e1dotq2[3],e1dotp[3],e2dotq1[3],e2dotp[3];
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) {
      e1dotq2[ix]=ZERO;
      e1dotp[ix] =ZERO;
      e2dotq1[ix]=ZERO;
      e2dotp[ix] =ZERO;
    }
    else {
      e1dotq2[ix] =_vectors[0][ix]*decay[2]->momentum();
      e1dotp[ix]  =_vectors[0][ix]*inpart.momentum();
      e2dotq1[ix] =_vectors[1][ix]*decay[1]->momentum();
      e2dotp[ix]  =_vectors[1][ix]*inpart.momentum();
    }
  }
  // the momentum dependent pieces of the matrix element
  Complex ii(0.,1.);
  Energy2 mpi2(sqr(decay[0]->mass())),meta2(sqr(inpart.mass())),
    mrho2(sqr(_rhomass)),
    t(mpi2+2.*((decay[0]->momentum())*(decay[1]->momentum()))),
    u(mpi2+2.*((decay[0]->momentum())*(decay[2]->momentum())));
  Energy q(sqrt(t)),pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
  complex<Energy2> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q=sqrt(u);pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy2> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-t-tgamma)),prop2(1./(mrho2-u-ugamma));
  complex<InvEnergy2> Dfact(_dconst[imode()]*(prop1*(pdotq2-meta2)
					      +prop2*(pdotq1-meta2)));
  complex<InvEnergy4> Efact(_econst[imode()]*(prop1+prop2));
  Complex e1dote2;
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      if(ix==1||iy==1) (*ME())(0,0,ix,iy)=0.;
      else {
	e1dote2=_vectors[0][ix].dot(_vectors[1][iy]);
	(*ME())(0,0,ix,iy) = 
	  Dfact*complex<Energy2>(e1dote2*q1dotq2-
				 e1dotq2[ix]*e2dotq1[iy])
	  -Efact*complex<Energy4>(-e1dote2*pdotq1*pdotq2
				  -e1dotp[ix]*e2dotp[iy]*q1dotq2
				  +e1dotq2[ix]*e2dotp[iy]*pdotq1
				  +e1dotp[ix]*e2dotq1[iy]*pdotq2);
      }
    }
  }
  /*
  double me(ME()->contract(rhoin).real());
  Energy M(inpart.mass()),M2(M*M);
  Energy2 s1(2.*(decay[1]->momentum()*decay[2]->momentum()));
  Energy2 s2(M2-2.*(inpart.momentum()*decay[1]->momentum()));
  Energy2 s3(M2-2.*(inpart.momentum()*decay[2]->momentum()));
  cout << "testing the matrix element " << (
   2*(2*(Dfact*conj(Dfact)).real() + 2*(Dfact*conj(Efact)).real()*M2 
      + (Efact*conj(Efact)).real()*M2*M2)*
      s1*s1 - 2*(Efact*conj(Efact)).real()*M2*s1*(M2 - s2)*
   (M2 - s3) +(Efact*conj(Efact)).real()*(M2 - s2)*(M2 - s2)*
   (M2-s3)*(M2-s3))/8. - me << endl;
  return me;
  */
  return ME()->contract(_rho).real();
}
 
double EtaPiGammaGammaDecayer::
threeBodyMatrixElement(const int imodeb, const Energy2 q2,const  Energy2 s3,
		       const Energy2 s2,const Energy2 s1,const Energy ,
		       const Energy ,const Energy ) const {
  // compute the prefactors
  Energy2 mrho2 = sqr(_rhomass);
  Energy q = sqrt(s3);
  Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  Complex ii(0.,1.);
  complex<Energy2> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q = sqrt(s2);
  pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy2> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-s3-tgamma)), prop2(1./(mrho2-s2-ugamma));
  Energy2 pdotq2(0.5*(q2-s3)), pdotq1(0.5*(q2-s2));
  complex<InvEnergy2> Dfact(_dconst[imodeb]*(prop1*(pdotq2-q2)+prop2*(pdotq1-q2)));
  complex<InvEnergy4> Efact(_econst[imodeb]*(prop1+prop2));
  InvEnergy4 D2 = (Dfact*conj(Dfact)).real();
  InvEnergy8 E2((Efact*conj(Efact)).real());
  InvEnergy6 ED((Efact*conj(Dfact)).real());
  return (2 * (2*D2 + 2*ED*q2 + E2*sqr(q2)) * sqr(s1)
	  - double(2*E2*q2*s1*(q2-s2)*(q2-s3))
	  + double(E2*sqr(q2-s2)*sqr(q2-s3))
	  )/8.;
}

WidthCalculatorBasePtr 
EtaPiGammaGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights; inweights.push_back(0.5); inweights.push_back(0.5);
  Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
  Energy wrho(getParticleData(ParticleID::rhoplus)->width());
  vector<Energy> inmass;  inmass.push_back(mrho);  inmass.push_back(mrho);
  vector<Energy> inwidth; inwidth.push_back(wrho); inwidth.push_back(wrho);
  vector<int> intype; intype.push_back(1); intype.push_back(2);
  vector<double> inpow(2,0.0);
  return new_ptr(ThreeBodyAllOnCalculator<EtaPiGammaGammaDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,
		  imode,_mpi,ZERO,ZERO));
}

void EtaPiGammaGammaDecayer::dataBaseOutput(ofstream & output, 
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":grhoomega " << _grhoomega*GeV << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV  << "\n";
  output << "newdef " << name() << ":grho " << _grho << "\n";
  output << "newdef " << name() << ":RhoMass " << _rhomass/MeV << "\n";
  output << "newdef " << name() << ":RhoWidth " << _rhowidth/MeV << "\n";
  output << "newdef " << name() << ":RatioFpiF8 " << _ratiofpif8 << "\n";
  output << "newdef " << name() << ":RatioFpiF0 " << _ratiofpif0 << "\n";
  output << "newdef " << name() << ":Theta " << _theta  << "\n";
  output << "newdef " << name() << ":EtaMax " << _etamax << "\n";
  output << "newdef " << name() << ":EtaPrimeMax " << _etapmax << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
