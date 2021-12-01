// -*- C++ -*-
//
// VectorMeson2MesonDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2MesonDecayer class.
//

#include "VectorMeson2MesonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig; 
using namespace ThePEG::Helicity;
 
void VectorMeson2MesonDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxweight_[ix]= mode(ix)->maxWeight();
  }
}

void VectorMeson2MesonDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxweight_.size() ||
     isize!=coupling_.size()) {
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMeson2MesonDecayer" << Exception::runerror;
  }
  // set up the integration channelsx
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in =  getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

VectorMeson2MesonDecayer::VectorMeson2MesonDecayer() {
  // don't generate intermediates
  generateIntermediates(false);
}
 
int VectorMeson2MesonDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(idbar==incoming_[ix]&&imode<0) {
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
  
void VectorMeson2MesonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxweight_ << coupling_;
}
  
void VectorMeson2MesonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxweight_ >> coupling_;
}
  
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMeson2MesonDecayer,DecayIntegrator>
describeHerwigVectorMeson2MesonDecayer("Herwig::VectorMeson2MesonDecayer", "HwVMDecay.so");

void VectorMeson2MesonDecayer::Init() {
  
  static ClassDocumentation<VectorMeson2MesonDecayer> documentation
    ("The VectorMeson2MesonDecayer class is designed to implement "
     "the decay of vector mesons to 2 scalar mesons via a current which is the "
     "difference of the momenta of the two scalars. The order of the scalar meson "
     "momenta does not matter as it only changes the sign of the matrix element.");

  static Command<VectorMeson2MesonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, pseudoscalars, coupling and max weight for a decay",
     &VectorMeson2MesonDecayer::setUpDecayMode, false);

  static Deleted<VectorMeson2MesonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");
  
}
 
void VectorMeson2MesonDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}
 
double VectorMeson2MesonDecayer::me2(const int,const Particle & part,
				     const tPDVector & ,
				     const vector<Lorentz5Momentum> & momenta,
				     MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  // difference of the momenta
  Lorentz5Vector<double> pdiff = (momenta[0]-momenta[1]) 
    * coupling_[imode()]/part.mass();
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) 
    (*ME())(ix,0,0)=vectors_[ix].dot(pdiff);
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = 4.*sqr(coupling_[imode()]*pcm/part.mass())/3.;
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}
 
bool VectorMeson2MesonDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					     double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()) idbar=dm.parent()->CC()->id();
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()) id1bar=(**pit).CC()->id();
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()) id2bar=(**pit).CC()->id();
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
  coupling=coupling_[imode];
  mecode=0;
  return order;
}

// output the setup information for the particle database
void VectorMeson2MesonDecayer::dataBaseOutput(ofstream & output,
					      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string VectorMeson2MesonDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
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
  coupling_.push_back(g);
  maxweight_.push_back(wgt);
  // success
  return "";
}
