// -*- C++ -*-
//
// PScalarLeptonNeutrinoDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarLeptonNeutrinoDecayer class.
//

#include "PScalarLeptonNeutrinoDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarLeptonNeutrinoDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  unsigned int iz(0),ix,iy;
  if(initialize()) {
    for(ix=0;ix<incoming_.size();++ix) {
      for(iy=0;iy<leptons_[ix];++iy) {
	maxWeight_[ix][iy] = mode(iz)->maxWeight();
	++iz;
      }
    }
  }
}

void PScalarLeptonNeutrinoDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=decayConstant_.size()||isize!=leptons_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in PScalarLeptonNeutrinoDecayer"
			  << Exception::abortnow;
  // create the integration channels
  tPDPtr nu[3]={getParticleData(ParticleID::nu_e),
		getParticleData(ParticleID::nu_mu),
		getParticleData(ParticleID::nu_tau)};
  tPDPtr nubar[3]={getParticleData(ParticleID::nu_ebar),
		   getParticleData(ParticleID::nu_mubar),
		   getParticleData(ParticleID::nu_taubar)};
  tPDPtr lep[3]={getParticleData(ParticleID::eminus),
		 getParticleData(ParticleID::muminus),
		 getParticleData(ParticleID::tauminus)};
  tPDPtr lepbar[3]={getParticleData(ParticleID::eplus),
		    getParticleData(ParticleID::muplus),
		    getParticleData(ParticleID::tauplus)};
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    for(unsigned int iy=0;iy<leptons_[ix];++iy) {
      tPDVector out(2);
      if(in->iCharge()>0) {
	out[0]=lepbar[iy];
	out[1]=nu[iy];
      }
      else {
	out[0]=lep[iy];
	out[1]=nubar[iy];
      }
      double wgt = maxWeight_[ix][iy];
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgt));
      addMode(mode);
    }
  }
}

int PScalarLeptonNeutrinoDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id()),id0bar(id0);
  if(parent->CC()){id0bar=-id0;}
  tPDVector::const_iterator pit = children.begin();
  int id;
  unsigned ilep(4);
  for(;pit!=children.end();++pit) {
    id=abs((**pit).id());
    if(id>=11&&id<=16&&id%2==0) ilep=(id-10)/2;
  }
  // find the channel we need
  bool found(false);
  int ichan(-1);
  unsigned int ix(0);
  do {
    if(id0   ==incoming_[ix]||id0bar==incoming_[ix]) {
      found=true;ichan+=ilep;
      cc=id0bar==incoming_[ix];
    }
    else {
      ichan+=leptons_[ix];
    }
    ++ix;
  }
  while (!found&&ix<incoming_.size());
  if(found) imode=ichan;
  return imode;
}


void PScalarLeptonNeutrinoDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << ounit(decayConstant_,GeV) << leptons_ << maxWeight_;
}

void PScalarLeptonNeutrinoDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> iunit(decayConstant_,GeV) >> leptons_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarLeptonNeutrinoDecayer,DecayIntegrator>
describeHerwigPScalarLeptonNeutrinoDecayer("Herwig::PScalarLeptonNeutrinoDecayer", "HwSMDecay.so");

void PScalarLeptonNeutrinoDecayer::Init() {

  static ClassDocumentation<PScalarLeptonNeutrinoDecayer> documentation
    ("The PScalarLeptonNeutrinoDecayer class is the base class for"
     " the implementation of leptonic decay of a pseudoscalar meson to a lepton"
     "and a neutrino.");

  static Command<PScalarLeptonNeutrinoDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, decay constant no of leptonic decays, max weights",
     &PScalarLeptonNeutrinoDecayer::setUpDecayMode, false);

  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceLeptons
    ("Leptons","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightElectron
    ("MaxWeightElectron","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightMuon
    ("MaxWeightMuon","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightTau
    ("MaxWeightTau","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceDecayConstant
    ("DecayConstant","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarLeptonNeutrinoDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(0);
  for(unsigned ix=0;ix<decay.size();++ix) {
    long id=decay[ix]->id();
    if(id<=-11&&id>=-16)    ianti=ix;
    else if(id>=11&&id<=16) iferm=ix;
  }
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double PScalarLeptonNeutrinoDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // work out which decay constant to use
  int icoup(0),id(abs(part.id()));
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(id==abs(incoming_[ix])) icoup=ix;
  }
  // find the particles
  unsigned int iferm(0),ianti(0);
  for(unsigned ix=0;ix<outgoing.size();++ix) {
    id=outgoing[ix]->id();
    if(id<=-11&&id>=-16)    ianti=ix;
    else if(id>=11&&id<=16) iferm=ix;
  }
  int idferm = outgoing[iferm]->id();
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // the prefactor
  Energy premass =  idferm%2==0 ? momenta[ianti].mass() : momenta[iferm].mass();
  InvEnergy pre = premass * 2.*decayConstant_[icoup]*SM().fermiConstant()/part.mass();
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(ispin[ianti+1]=0;ispin[ianti+1]<2;++ispin[ianti+1]) {
    for(ispin[iferm+1]=0;ispin[iferm+1]<2;++ispin[iferm+1]) {
      (*ME())(ispin)= idferm%2==0 ? 
	Complex(pre*wave_[ispin[ianti+1]].rightScalar(wavebar_[ispin[iferm+1]])) :
	Complex(pre*wave_[ispin[ianti+1]].leftScalar( wavebar_[ispin[iferm+1]])) ;
    }
  }
  double me = 0.5*ME()->contract(rho_).real();
  // test of the matrix element
  // Energy mass = idferm%2==0 ? momenta[ianti].mass() : momenta[iferm].mass();
  // double test = 0.5*sqr(decayConstant_[icoup]*SM().fermiConstant()*2.*mass/part.mass())*
  //   (sqr(part.mass())-sqr(mass));
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << (me-test)/(me+test) << endl;
  return me;
}

void PScalarLeptonNeutrinoDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << decayConstant_[ix]/MeV << " " << leptons_[ix];
    for(unsigned int iy=0;iy<3;++iy) output << " " << maxWeight_[ix][iy];
    output << "\n";
  }
  if(header) 
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalarLeptonNeutrinoDecayer::setUpDecayMode(string arg) {
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
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // no of leptons
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int lepton = stoi(stype);
  // and the maximum weights
  array<double,3> wgt;
  for(unsigned int iy=0;iy<3;++iy) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    wgt[iy] = stof(stype);
  }
  // store the information
  incoming_.push_back(in);
  decayConstant_.push_back(g*MeV);
  leptons_.push_back(lepton);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
