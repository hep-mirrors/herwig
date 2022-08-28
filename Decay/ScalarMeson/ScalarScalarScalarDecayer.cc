// -*- C++ -*-
//
// ScalarScalarScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarScalarScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix) ) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void ScalarScalarScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in ScalarScalarScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int ScalarScalarScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  // loop over the modes and see if this is one of them
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==incoming_[ix]&&imode<0) {
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

void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,MeV) << incoming_ << outgoing_ << maxWeight_;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,MeV) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarScalarScalarDecayer,DecayIntegrator>
describeHerwigScalarScalarScalarDecayer("Herwig::ScalarScalarScalarDecayer", "HwSMDecay.so");

void ScalarScalarScalarDecayer::Init() {

  static ClassDocumentation<ScalarScalarScalarDecayer> documentation
    ("The ScalarScalarScalarDecayer class is designed for the"
     " decay of a scalar meson to two scalar mesons including off-shell effects");

  static Command<ScalarScalarScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing scalars, coupling(MeV) and max weight for a decay",
     &ScalarScalarScalarDecayer::setUpDecayMode, false);
  
  static Deleted<ScalarScalarScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

}

void ScalarScalarScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double ScalarScalarScalarDecayer::me2(const int,const Particle & part,
				      const tPDVector &,
				      const vector<Lorentz5Momentum> &,
				      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  double fact(coupling_[imode()]/part.mass());
  (*ME())(0,0,0) = fact;
  return sqr(fact);
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0);
  bool order(true);
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
  coupling=coupling_[imode]/dm.parent()->mass();
  itype = 6;
  return order;
}

// output the setup information for the particle database
void ScalarScalarScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]/MeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string ScalarScalarScalarDecayer::setUpDecayMode(string arg) {
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
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 0";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
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
  coupling_.push_back(g*MeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
