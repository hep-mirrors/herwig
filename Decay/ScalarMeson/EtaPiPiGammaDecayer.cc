// -*- C++ -*-
//
// EtaPiPiGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiGammaDecayer class.
//

#include "EtaPiPiGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/Decay/FormFactors/OmnesFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


DescribeClass<EtaPiPiGammaDecayer,DecayIntegrator>
describeHerwigEtaPiPiGammaDecayer("Herwig::EtaPiPiGammaDecayer",
				  "HwSMDecay.so");
HERWIG_INTERPOLATOR_CLASSDESC(EtaPiPiGammaDecayer,double,Energy)



void EtaPiPiGammaDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_maxweight.size();++ix)
      _maxweight[ix]=mode(ix)->maxWeight();
  }
}

EtaPiPiGammaDecayer::EtaPiPiGammaDecayer() 
  : _incoming(2), _coupling(2), _maxweight(2), _option(2) {
  // the pion decay constant
  _fpi=130.7*MeV;
  // the rho mass
  _mrho=0.7711*GeV;
  _rhowidth=0.1492*GeV;
  // the constants for the omnes function form
  _aconst=0.5/_mrho/_mrho;
  _cconst=1.0;
  // use local values of the parameters
  _localparameters=true;
  // the modes
  // eta decay
  _incoming[0] = 221; 
  _option[0] = 3; 
  _coupling[0] = 5.060e-3; 
  _maxweight[0] = 3.95072; 
  // eta' decay
  _incoming[1] = 331; 
  _option[1] = 3; 
  _coupling[1] = 4.278e-3; 
  _maxweight[1] = 3.53141; 
  _rhoconst=0.;
  _mpi=ZERO;
  // intermediates
  generateIntermediates(false);
}

void EtaPiPiGammaDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()||isize!=_option.size()||isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "EtaPiPiGammaDecayer::doinit()" << Exception::abortnow;
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters) {
    _mrho=rho->mass();
    _rhowidth=rho->width();
  }
  _mpi=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(_mrho,_mpi,_mpi));
  _rhoconst=sqr(_mrho)*_rhowidth/pow<3,1>(pcm);
  // set up the modes
  tPDVector out = {getParticleData(ParticleID::piplus),
		   getParticleData(ParticleID::piminus),
		   getParticleData(ParticleID::gamma)};
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_coupling.size();++ix) {
    tPDPtr in = getParticleData(_incoming[ix]);
    mode = new_ptr(PhaseSpaceMode(in,out,_maxweight[ix]));
    mode->addChannel((PhaseSpaceChannel(mode),0,rho,0,3,1,1,1,2));
    addMode(mode);
  }
}

int EtaPiPiGammaDecayer::modeNumber(bool & cc,tcPDPtr parent,
				    const tPDVector & children) const {
  int imode(-1);
  // check number of external particles
  if(children.size()!=3){return imode;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),ngamma(0);
  tPDVector::const_iterator pit = children.begin();
  int id;
  for(;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::gamma)   ++ngamma;
  }
  if(!(npip==1&&npim==1&&ngamma==1)) return imode;
  unsigned int ix(0);
  id=parent->id();
  do{if(id==_incoming[ix]){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void EtaPiPiGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_fpi,MeV) << _incoming << _coupling << _maxweight << _option 
     << ounit(_aconst,1/MeV2) << _cconst <<ounit(_mrho,MeV) << ounit(_rhowidth,MeV) 
     << _rhoconst << ounit(_mpi,MeV) << _localparameters << omnesFunction_;
}

void EtaPiPiGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_fpi,MeV) >> _incoming >> _coupling >> _maxweight >> _option 
     >> iunit(_aconst,1/MeV2) >> _cconst >>iunit(_mrho,MeV) >> iunit(_rhowidth,MeV) 
     >> _rhoconst >> iunit(_mpi,MeV) >> _localparameters >> omnesFunction_;
}

void EtaPiPiGammaDecayer::Init() {

  static ClassDocumentation<EtaPiPiGammaDecayer> documentation
    ("The EtaPiPiGammaDecayer class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma",
     "The decays of $\\eta,\\eta'\\to\\pi^+\\pi^-\\gamma$ were simulated"
     " using the matrix elements from \\cite{Venugopal:1998fq,Holstein:2001bt}",
     "\\bibitem{Venugopal:1998fq} E.~P.~Venugopal and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 57} (1998) 4397 [arXiv:hep-ph/9710382].\n"
     "%%CITATION = PHRVA,D57,4397;%%\n"
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Reference<EtaPiPiGammaDecayer,OmnesFunction> interfaceOmnesFunction
    ("OmnesFunction",
     "Omnes function for the matrix element",
     &EtaPiPiGammaDecayer::omnesFunction_, false, false, true, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiGammaDecayer::_fpi, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static ParVector<EtaPiPiGammaDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &EtaPiPiGammaDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &EtaPiPiGammaDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &EtaPiPiGammaDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiGammaDecayer::_mrho, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Parameter<EtaPiPiGammaDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_cconst, 1.0, -10., 10.,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_aconst, 1./GeV2, 0.8409082/GeV2, ZERO,
     10./GeV2,
     false, false, false);

  static ParVector<EtaPiPiGammaDecayer,int> interfaceOption
    ("Option",
     "The form of the prefactor 0 is a VMD model using M Gamma for the width term,"
     "1 is a VMD model using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiGammaDecayer::_option,
     0, 0, 0, 0, 4, false, false, true);
}
void EtaPiPiGammaDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  VectorWaveFunction::constructSpinInfo(_vectors,decay[2],
					outgoing,true,true);
}

double EtaPiPiGammaDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  _vectors.resize(3);
  for(unsigned int ix=0;ix<3;ix+=2) {
    _vectors[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // prefactor for the matrix element
  complex<InvEnergy3> pre(_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi));
  Lorentz5Momentum ppipi(momenta[0]+momenta[1]);
  ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imode()]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*_mrho*pcm*pcm*pcm*_rhoconst/q2));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imode()]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(_option[imode()]==2 || _option[imode()]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*q2)/omnesFunction_->D(q2));
  }
  pre = pre*fact;
  LorentzPolarizationVector epstemp(pre*Helicity::epsilon(momenta[0],
							  momenta[1],
							  momenta[2]));
  // compute the matrix element
  vector<unsigned int> ispin(4,0);
  for(ispin[3]=0;ispin[3]<3;++ispin[3]) {
    if(ispin[3]==1) (*ME())(ispin)=0.;
    else            (*ME())(ispin)=epstemp.dot(_vectors[ispin[3]]);
  }
  // contract the whole thing
  return ME()->contract(_rho).real();
}

double EtaPiPiGammaDecayer::
threeBodyMatrixElement(const int imodeb,const Energy2 ,const  Energy2 s3,const 
		       Energy2 s2,const Energy2 s1,const Energy ,
		       const Energy ,const Energy ) const {
  complex<InvEnergy3> pre(_coupling[imodeb]*2.*sqrt(2.)/pow<3,1>(_fpi));
  Energy q(sqrt(s3));
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imodeb]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*_mrho*pcm*pcm*pcm*_rhoconst/s3));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imodeb]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(_option[imodeb]==2||_option[imodeb]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*s3)/omnesFunction_->D(s3));
  }
  pre =pre*fact;
  InvEnergy6 factor((pre*conj(pre)).real());
  Energy2 mpi2(_mpi*_mpi);
  return factor*((-mpi2*(-2*mpi2+s1+s2)*(-2*mpi2+s1+s2)+(mpi2-s1)*(mpi2-s2)*s3)/4.);
}

WidthCalculatorBasePtr 
EtaPiPiGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights(1,1.);
  vector<Energy> inmass(1,getParticleData(ParticleID::rho0)->mass());
  vector<Energy> inwidth(1,getParticleData(ParticleID::rho0)->width());
  vector<int> intype(1,1);
  vector<double> inpow(1,0.0);
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator<EtaPiPiGammaDecayer>
		   (inweights,intype,inmass,inwidth,inpow,*this,imode,_mpi,_mpi,ZERO)));
  return output;
}

void EtaPiPiGammaDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":fpi             " << _fpi/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << _mrho/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << _rhowidth/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":OmnesC          " << _cconst          << "\n";
  output << "newdef " << name() << ":OmnesA          " << _aconst*GeV2     << "\n";
  for(unsigned int ix=0;ix<2;++ix) {
    output << "newdef " << name() << ":Incoming    " << ix << "  " 
	   << _incoming[ix]    << "\n";
    output << "newdef " << name() << ":Coupling    " << ix << "  " 
	   << _coupling[ix]    << "\n";
    output << "newdef " << name() << ":MaxWeight   " << ix << "  " 
	   << _maxweight[ix]   << "\n";
    output << "newdef " << name() << ":Option      " << ix << "  " 
	   << _option[ix]      << "\n";
  }
  omnesFunction_->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":OmnesFunction " << omnesFunction_->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
