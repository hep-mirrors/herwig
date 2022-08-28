// -*- C++ -*-
//
// VectorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorPScalarDecayer class.
//

#include "VectorMesonVectorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()  || 
     isize!=maxweight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorPScalarDecayer" << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

VectorMesonVectorPScalarDecayer::VectorMesonVectorPScalarDecayer()  {
  // intermediates
  generateIntermediates(false);
}

int VectorMesonVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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


void VectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxweight_ << ounit(coupling_,1/MeV);
}

void VectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxweight_ >> iunit(coupling_,1/MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonVectorPScalarDecayer,DecayIntegrator>
describeHerwigVectorMesonVectorPScalarDecayer("Herwig::VectorMesonVectorPScalarDecayer", "HwVMDecay.so");

void VectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorPScalarDecayer> documentation
    ("The VectorMesonVectorPScalarDecayer class is designed for the "
     "decay of a vector meson to another vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static Command<VectorMesonVectorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, pesudovector), coupling(1/MeV) and max weight for a decay",
     &VectorMesonVectorPScalarDecayer::setUpDecayMode, false);

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceOutcomingVector
    ("OutgoingVector","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceOutcomingPScalar
    ("OutgoingPScalar","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMesonVectorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[0],
					outgoing,true,decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double VectorMesonVectorPScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  // is the vector massless
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  
  // compute the matrix element
  LorentzPolarizationVector vtemp;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=0.;
    }
    else {
      vtemp=coupling_[imode()]/part.mass()*
	epsilon(part.momentum(),vectors_[1][ix],momenta[0]);
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=vectors_[0][iy].dot(vtemp);
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = sqr(coupling_[imode()]*pcm)*2./3.;
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
    if(id==incoming_[ix]) {
      if(id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second) {
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
  coupling = coupling_[imode]*dm.parent()->mass();  
  mecode = 1;
  return order;
}

// output the setup info for the particle database
void VectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*MeV << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMesonVectorPScalarDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy g = stof(stype)/MeV;
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
