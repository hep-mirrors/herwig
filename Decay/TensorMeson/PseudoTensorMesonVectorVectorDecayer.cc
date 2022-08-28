// -*- C++ -*-
//
// PseudoTensorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoTensorMesonVectorVectorDecayer class.
//

#include "PseudoTensorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoTensorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void PseudoTensorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters PseudoTensorMesonVectorVectorDecayer" 
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

int PseudoTensorMesonVectorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
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
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void PseudoTensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1./GeV);
}

void PseudoTensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoTensorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigPseudoTensorMesonVectorVectorDecayer("Herwig::PseudoTensorMesonVectorVectorDecayer", "HwTMDecay.so");

void PseudoTensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<PseudoTensorMesonVectorVectorDecayer> documentation
    ("The PseudoTensorMesonVectorVectorDecayer class implements the"
     " decay of a pseudotensor meson to two vector mesons");

  static Command<PseudoTensorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &PseudoTensorMesonVectorVectorDecayer::setUpDecayMode, false);

}
void PseudoTensorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,
					  decay[ix]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double PseudoTensorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  // check for photons
  bool photon[2] = {outgoing_[imode()].first ==ParticleID::gamma,
		    outgoing_[imode()].second==ParticleID::gamma};
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  Energy2 denom[2];
  for(unsigned int iy=0;iy<2;++iy) {
    denom[iy] = momenta[iy]*part.momentum()-momenta[iy].mass()*part.mass();
    vectors_[iy].resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(photon[iy] && ix==1) continue;
      vectors_[iy][ix] = HelicityFunctions::polarizationVector(-momenta[iy],ix,Helicity::outgoing);
    }
  }
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass());
  InvEnergy2 c1 = 0.5*(sqr(part.mass())+sqr(momenta[0].mass())-sqr(momenta[1].mass()))/sqr(part.mass()*pcm);
  InvEnergy2 c2 = 0.5*(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass()))/sqr(part.mass()*pcm);
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzPolarizationVectorE v0 = tensors_[inhel].postDot(momenta[0]-momenta[1]);
    LorentzVector<complex<Energy3> > v1 = epsilon(v0,momenta[0],momenta[1]);
    for(unsigned int vhel1=0;vhel1<3;++vhel1) {
      LorentzVector<complex<Energy2> > v3 = epsilon(part.momentum(),v0,vectors_[0][vhel1]);
      complex<InvEnergy> d1 = c1*vectors_[0][vhel1]*momenta[1];
      complex<Energy3> d3 =  vectors_[0][vhel1]*v1;
      for(unsigned int vhel2=0;vhel2<3;++vhel2) {
	if ( (photon[0] && vhel1==1) || (photon[1] && vhel2==1))
	  (*ME())(inhel,vhel1,vhel2)=0.;
	else {
	  complex<InvEnergy> d2 = c2*vectors_[1][vhel2]*momenta[0];
	  (*ME())(inhel,vhel1,vhel2)=Complex(fact*(v3.dot(vectors_[1][vhel2])
						   - d1*v1.dot(vectors_[1][vhel2]) - d2*d3));
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test = 16./15.*sqr(coupling_[imode()]*pcm);
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool PseudoTensorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=9;
  return order;
}

void PseudoTensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
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

string PseudoTensorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
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
  maxWeight_.push_back(wgt);
  // success
  return "";
}
