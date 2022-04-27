// -*- C++ -*-
//
// TensorMeson2PScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMeson2PScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMeson2PScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMeson2PScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
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

int TensorMeson2PScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV);
}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMeson2PScalarDecayer,DecayIntegrator>
describeHerwigTensorMeson2PScalarDecayer("Herwig::TensorMeson2PScalarDecayer", "HwTMDecay.so");

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The TensorMeson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static Command<TensorMeson2PScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, scalars, coupling(1/GeV) and max weight for a decay",
     &TensorMeson2PScalarDecayer::setUpDecayMode, false);
  
  static Deleted<TensorMeson2PScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");
}

void TensorMeson2PScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// matrix elememt for the process
double TensorMeson2PScalarDecayer::me2(const int,const Particle & part,
				       const tPDVector &,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  // calculate the matrix element
  for(unsigned int ix=0;ix<5;++ix) {
    (*ME())(ix,0,0) = Complex(coupling_[imode()]/part.mass()*
			      ((tensors_[ix]*momenta[1])*momenta[0]));
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = Energy4(pow<4,1>(2*pcm))*sqr( coupling_[imode()]/part.mass())/120.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=7;
  return order;
}

void TensorMeson2PScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMeson2PScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
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
  coupling_.push_back(g/GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
