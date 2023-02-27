// -*- C++ -*-
//
// PScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorVectorDecayer class.
//
#include "PScalarVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix] = mode(ix)->maxWeight();
  }
}

void PScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  = getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  unsigned int ix(0);
  int imode(-1);
  do {
    if(incoming_[ix]==id0) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) {
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
  while(imode<0&&ix<incoming_.size());
  return imode;
}

void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/GeV) << incoming_ << outgoing_ << maxweight_;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/GeV) >> incoming_ >> outgoing_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarVectorVectorDecayer,DecayIntegrator>
describeHerwigPScalarVectorVectorDecayer("Herwig::PScalarVectorVectorDecayer", "HwSMDecay.so");

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The PScalarVectorVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static Command<PScalarVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &PScalarVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->id()==ParticleID::gamma;
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,photon[ix]);
}

double PScalarVectorVectorDecayer::me2(const int,const Particle & part,
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
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon[0] && ix==1) continue;
    for(unsigned int iy=0;iy<3;++iy) {
      if(photon[1] && iy==1) continue;
      (*ME())(0,ix,iy)=Complex(fact*epsilon(vectors_[0][ix],momenta[1],
					    vectors_[1][iy])*momenta[0]);
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 2.*sqr(fact*part.mass())*
  //   sqr(Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass()));
  // cerr << "testing the matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << (output-test)/(output+test) << "\n";
  return output;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  while(imode<0&&ix<incoming_.size());
  coupling=coupling_[imode]*dm.parent()->mass();
  itype = 3;
  return order;
}

// output the setup info for the particle database
void PScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxweight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PScalarVectorVectorDecayer::setUpDecayMode(string arg) {
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
  coupling_.push_back(g/GeV);
  maxweight_.push_back(wgt);
  // success
  return "";
}
