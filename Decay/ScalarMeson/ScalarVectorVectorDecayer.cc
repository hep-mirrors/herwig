// -*- C++ -*-
//
// ScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarVectorVectorDecayer class.
//

#include "ScalarVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void ScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "ScalarVectorVectorDecayerDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in     =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first ),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int ScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  cc = false;
  // check that at least some modes exist
  // must be two outgoing particles
  if(incoming_.size()==0||children.size()!=2) return -1;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  int ccid0(parent     ->CC() ? parent     ->CC()->id() : parent     ->id());
  int ccid1(children[0]->CC() ? children[0]->CC()->id() : children[0]->id());
  int ccid2(children[1]->CC() ? children[1]->CC()->id() : children[1]->id());
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  int imode(-1);
  do {
    if(incoming_[ix]==id0) {
      if((outgoing_[ix].first==id1&&outgoing_[ix].second==id2)||
	 (outgoing_[ix].first==id2&&outgoing_[ix].second==id1)) {
	imode=ix;
	break;
      }
    }
    if(incoming_[ix]==ccid0) {
      if((outgoing_[ix].first==ccid1&&outgoing_[ix].second==ccid2)||
	 (outgoing_[ix].first==ccid2&&outgoing_[ix].second==ccid1)) {
	imode=ix;
	cc=true;
	break;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

void ScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,GeV) << incoming_ << outgoing_ << maxWeight_;
}

void ScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,GeV) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarVectorVectorDecayer,DecayIntegrator>
describeHerwigScalarVectorVectorDecayer("Herwig::ScalarVectorVectorDecayer", "HwSMDecay.so");

void ScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<ScalarVectorVectorDecayer> documentation
    ("The ScalarVectorVectorDecayer class is designed for"
     " the decay of a scalar meson to two spin-1 particles.");

  static Command<ScalarVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &ScalarVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<ScalarVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");
}

void ScalarVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double ScalarVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],ihel,Helicity::outgoing);
    }
  }
  // now compute the matrix element
  double fact(coupling_[imode()]/part.mass());
  Energy2 p1p2(momenta[0]*momenta[1]);
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix) {
    if(photon[0] && ix==1) continue;
    for(iy=0;iy<3;++iy) {
      if(photon[1] && iy==1) continue;
      (*ME())(0,ix,iy)=Complex(fact*(vectors_[0][ix].dot(vectors_[1][iy])-
				     (vectors_[1][iy]*momenta[0])*
				     (vectors_[0][ix]*momenta[1])/(p1p2-momenta[0].mass()*momenta[1].mass())));
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = sqr(coupling_[imode()]/part.mass())*
  //   (2.*sqr(pcm*part.mass())+3.*sqr(momenta[0].mass()*momenta[1].mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  return me;
}

// output the setup info for the particle database
void ScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					      double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do {
    if(incoming_[ix]==id) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  coupling=coupling_[imode]/dm.parent()->mass();
  itype = 12;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

string ScalarVectorVectorDecayer::setUpDecayMode(string arg) {
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
  coupling_.push_back(g*GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
