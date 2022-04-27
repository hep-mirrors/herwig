// -*- C++ -*-
//
// PScalarPScalarVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarPScalarVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void PScalarPScalarVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in PScalarPScalarVectorDecayer" 
  			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PScalarPScalarVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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


void PScalarPScalarVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxWeight_;
}

void PScalarPScalarVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarPScalarVectorDecayer,DecayIntegrator>
describeHerwigPScalarPScalarVectorDecayer("Herwig::PScalarPScalarVectorDecayer", "HwSMDecay.so");

void PScalarPScalarVectorDecayer::Init() {

  static ClassDocumentation<PScalarPScalarVectorDecayer> documentation
    ("The PScalarPScalarVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");
  static Command<PScalarPScalarVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling (dimensionless) and max weight for a decay",
     &PScalarPScalarVectorDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarPScalarVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceOutcomingP
    ("OutgoingP","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceOutcomingV
    ("OutgoingV","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarPScalarVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],
					outgoing,true,false);
}

double PScalarPScalarVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  if(meopt==Initialize) {
     ScalarWaveFunction::
       calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  vectors_.resize(3);
  bool massless = outgoing[1]->id()==ParticleID::gamma;
  for(unsigned int ix=0;ix<3;++ix) {
    if(massless && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Lorentz5Momentum psum(part.momentum()+momenta[0]);
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(0,0,ix) = Complex(coupling_[imode()]/part.mass()*(vectors_[ix]*psum));
  }
  double me=ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 4.*sqr(coupling_[imode()]*pcm/momenta[1].mass());
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarPScalarVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
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
  unsigned int ix(0); bool order(true);
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
  coupling=coupling_[imode];
  mecode=10;
  return order;
}

// output the setup information for the particle database
void PScalarPScalarVectorDecayer::dataBaseOutput(ofstream & output,
						 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix] << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PScalarPScalarVectorDecayer::setUpDecayMode(string arg) {
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
