// -*- C++ -*-
//
// EtaPiPiFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiFermionsDecayer class.
//

#include "EtaPiPiFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/Decay/FormFactors/OmnesFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


DescribeClass<EtaPiPiFermionsDecayer,DecayIntegrator>
describeHerwigEtaPiPiFermionsDecayer("Herwig::EtaPiPiFermionsDecayer",
				  "HwSMDecay.so");
HERWIG_INTERPOLATOR_CLASSDESC(EtaPiPiFermionsDecayer,double,Energy)



void EtaPiPiFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
      maxWeight_[ix]=mode(ix)->maxWeight();
    }
  }
}

EtaPiPiFermionsDecayer::EtaPiPiFermionsDecayer() 
  : fPi_(130.7*MeV), mRho_(0.7711*GeV), rhoWidth_(0.1492*GeV) {
  // the constants for the omnes function form
  aConst_=0.5/mRho_/mRho_;
  cConst_=1.0;
  // use local values of the parameters
  localParameters_=true;
  // the modes
  mPi_=ZERO;
  // intermediates
  generateIntermediates(false);
}

void EtaPiPiFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=coupling_.size()||isize!=option_.size()||isize!=lepton_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "EtaPiPiFermionsDecayer::doinit()" << Exception::abortnow;
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!localParameters_) {
    mRho_=rho->mass();
    rhoWidth_=rho->width();
  }
  mPi_=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(mRho_,mPi_,mPi_));
  rhoConst_=sqr(mRho_)*rhoWidth_/pow<3,1>(pcm);
  // set up the modes
  tPDPtr gamma = getParticleData(ParticleID::gamma);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<coupling_.size();++ix) {
    tPDVector out = {getParticleData(ParticleID::piplus),
		     getParticleData(ParticleID::piminus),
		     getParticleData( lepton_[ix]),
		     getParticleData(-lepton_[ix])};
    tPDPtr in = getParticleData(incoming_[ix]);
    mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,rho,0,gamma,1,1,1,2,2,3,2,4));
    newChannel.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
    addMode(mode);
  }
}

int EtaPiPiFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
				    const tPDVector & children) const {
  int imode(-1);
  // check number of external particles
  if(children.size()!=4){return imode;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),nl(0);
  int il(0);
  tPDVector::const_iterator pit = children.begin();
  int id;
  for(;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else {
      il = abs(id);
      ++nl;
    }
  }
  if(!(npip==1&&npim==1&&nl==2)) return imode;
  unsigned int ix(0);
  id=parent->id();
  do {
    if(id==incoming_[ix] && il==lepton_[ix]) imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void EtaPiPiFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << incoming_ << coupling_ << maxWeight_ << lepton_ << option_ 
     << ounit(aConst_,1/MeV2) << cConst_ <<ounit(mRho_,MeV) << ounit(rhoWidth_,MeV) 
     << rhoConst_ << ounit(mPi_,MeV) << localParameters_ << omnesFunction_;
}

void EtaPiPiFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> incoming_ >> coupling_ >> maxWeight_ >> lepton_ >> option_ 
     >> iunit(aConst_,1/MeV2) >> cConst_ >>iunit(mRho_,MeV) >> iunit(rhoWidth_,MeV) 
     >> rhoConst_ >> iunit(mPi_,MeV) >> localParameters_ >> omnesFunction_;
}

void EtaPiPiFermionsDecayer::Init() {

  static ClassDocumentation<EtaPiPiFermionsDecayer> documentation
    ("The EtaPiPiFermionsDecayer class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma",
     "The decays of $\\eta,\\eta'\\to\\pi^+\\pi^-\\gamma$ were simulated"
     " using the matrix elements from \\cite{Venugopal:1998fq,Holstein:2001bt}",
     "\\bibitem{Venugopal:1998fq} E.~P.~Venugopal and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 57} (1998) 4397 [arXiv:hep-ph/9710382].\n"
     "%%CITATION = PHRVA,D57,4397;%%\n"
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Reference<EtaPiPiFermionsDecayer,OmnesFunction> interfaceOmnesFunction
    ("OmnesFunction",
     "Omnes function for the matrix element",
     &EtaPiPiFermionsDecayer::omnesFunction_, false, false, true, false, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiFermionsDecayer::fPi_, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static Command<EtaPiPiFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, lepton, option for VMD, coupling and maxweight."
     "There are three options for for VMD, 0 is a VMD model using M Gamma for the width term,"
     " 1 is a VMD model using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiFermionsDecayer::setUpDecayMode, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiFermionsDecayer::mRho_, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiFermionsDecayer::rhoWidth_, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<EtaPiPiFermionsDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiFermionsDecayer::localParameters_, true, false, false);
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

  static Parameter<EtaPiPiFermionsDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiFermionsDecayer::cConst_, 1.0, -10., 10.,
     false, false, false);

  static Parameter<EtaPiPiFermionsDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiFermionsDecayer::aConst_, 1./GeV2, 0.8409082/GeV2, ZERO,
     10./GeV2,
     false, false, false);

}

void EtaPiPiFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[2],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[3],outgoing,true);
}

double EtaPiPiFermionsDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,
					 PDT::Spin1Half,PDT::Spin1Half)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  wavebar_.resize(2);
  wave_   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[2],ihel,Helicity::outgoing);
    wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[3],ihel,Helicity::outgoing);
  }
  Lorentz5Momentum pff(momenta[2]+momenta[3]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // prefactor for the matrix element
  complex<InvEnergy4> pre(part.mass()*coupling_[imode()]*2.*sqrt(2.)/(fPi_*fPi_*fPi_)*
			  sqrt(SM().alphaEM(mff2)*4.*Constants::pi)/mff2);
  Lorentz5Momentum ppipi(momenta[0]+momenta[1]);
  ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(option_[imode()]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,mPi_,mPi_));
    Complex resfact(q2/(mRho_*mRho_-q2-ii*mRho_*pcm*pcm*pcm*rhoConst_/q2));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(option_[imode()]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,mPi_,mPi_));
    Complex resfact(q2/(mRho_*mRho_-q2-ii*pcm*pcm*pcm*rhoConst_/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(option_[imode()]==2 || option_[imode()]==3) {
    fact=(1.-cConst_+cConst_*(1.+aConst_*q2)/omnesFunction_->D(q2));
  }
  pre = pre*fact;
  LorentzVector<complex<InvEnergy> >
    epstemp(pre*Helicity::epsilon(momenta[0],momenta[1],pff));
  // compute the matrix element
  vector<unsigned int> ispin(5,0);
  for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
    for(ispin[4]=0;ispin[4]<2;++ispin[4]) {
      LorentzPolarizationVectorE fcurrent = wave_[ispin[4]].vectorCurrent(wavebar_[ispin[3]]);
      (*ME())(ispin) = Complex(epstemp.dot(fcurrent));
    }
  }
  // contract the whole thing
  return ME()->contract(rho_).real();
}

void EtaPiPiFermionsDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":fpi             " << fPi_/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << mRho_/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << rhoWidth_/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << localParameters_ << "\n";
  output << "newdef " << name() << ":OmnesC          " << cConst_          << "\n";
  output << "newdef " << name() << ":OmnesA          " << aConst_*GeV2     << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << lepton_[ix] << " " << option_[ix] << " "  
	   << coupling_[ix] << " " << maxWeight_[ix] << "\n";
  }
  omnesFunction_->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":OmnesFunction " << omnesFunction_->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string EtaPiPiFermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "Outgoing fermion with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Outgoing fermion with id " + std::to_string(out) + "does not have spin 1/2";
  // vmd option  
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDinclude = stoi(stype);
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  lepton_.push_back(out);
  coupling_.push_back(g);
  option_.push_back(VMDinclude);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
