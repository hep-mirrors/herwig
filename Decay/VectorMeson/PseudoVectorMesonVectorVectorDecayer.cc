// -*- C++ -*-
//
// PseudoVectorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoVectorMesonVectorVectorDecayer class.
//

#include "PseudoVectorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoVectorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix]=mode(ix)->maxWeight();
  }
}

void PseudoVectorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()  ||
     isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "PseudoVectorMesonVectorVectorDecayer" << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

PseudoVectorMesonVectorVectorDecayer::PseudoVectorMesonVectorVectorDecayer() {
  // intermediates
  generateIntermediates(false);
}

int PseudoVectorMesonVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
						 const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void PseudoVectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_;
}

void PseudoVectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoVectorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigPseudoVectorMesonVectorVectorDecayer("Herwig::PseudoVectorMesonVectorVectorDecayer", "HwVMDecay.so");

void PseudoVectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<PseudoVectorMesonVectorVectorDecayer> documentation
    ("The PseudoVectorMesonVectorVectorDecayer class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static Command<PseudoVectorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling and max weight for a decay",
     &PseudoVectorMesonVectorVectorDecayer::setUpDecayMode, false);

}

string PseudoVectorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
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
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void PseudoVectorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix+1],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double PseudoVectorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix) 
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix+1].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix+1][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],
								   ihel,Helicity::outgoing);
    }
  }
  // work out the dot products we need for the matrix element
  complex<InvEnergy> fact(coupling_[imode()]/part.mass());
  Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  					   momenta[1].mass());
  InvEnergy3 kinFact(1./part.mass()/sqr(pcm));
  Energy E1(0.5*(sqr(part.mass()) + sqr(momenta[0].mass()) - sqr(momenta[1].mass()))/part.mass());
  Energy E2(0.5*(sqr(part.mass()) - sqr(momenta[0].mass()) + sqr(momenta[1].mass()))/part.mass());
  for(unsigned int ih0=0;ih0<3;++ih0) {
    for(unsigned int ih1=0;ih1<3;++ih1) {
      if(photon[0] && ih1==1) continue;
      LorentzPolarizationVectorE       v0 = epsilon(part.momentum(),vectors_[0][ih0],vectors_[1][ih1]);
      LorentzVector<complex<Energy2> > v1 = epsilon(vectors_[0][ih0],momenta[0],momenta[1]);
      complex<Energy2> d1 = v1*vectors_[1][ih1];
      complex<Energy> d3 = momenta[1]*vectors_[1][ih1];
      for(unsigned int ih2=0;ih2<3;++ih2) {
	complex<Energy2> d2 = v1*vectors_[2][ih2];
	complex<Energy> d4 = momenta[0]*vectors_[2][ih2];
	if(photon[1] && ih2==1) continue;
	(*ME())(ih0,ih1,ih2)= Complex(fact*(v0*vectors_[2][ih2]
					    -kinFact*(E1*d2*d3+E2*d1*d4)));
	
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element;
  // double test = 2./3.*sqr(coupling_[imode()]);
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool PseudoVectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						    double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling = coupling_[imode]; 
  mecode = 5;
  return order; 
}

// output the setup information for the particle database
void PseudoVectorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
