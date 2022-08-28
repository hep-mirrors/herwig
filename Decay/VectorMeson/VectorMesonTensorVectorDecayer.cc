// -*- C++ -*-
//
// VectorMesonTensorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonTensorVectorDecayer class.
//

#include "VectorMesonTensorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonTensorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void VectorMesonTensorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters VectorMesonTensorVectorDecayer" 
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

int VectorMesonTensorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
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


void VectorMesonTensorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void VectorMesonTensorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonTensorVectorDecayer,DecayIntegrator>
describeHerwigVectorMesonTensorVectorDecayer("Herwig::VectorMesonTensorVectorDecayer", "HwTMDecay.so");

void VectorMesonTensorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonTensorVectorDecayer> documentation
    ("The VectorMesonTensorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static Command<VectorMesonTensorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(GeV) and max weight for a decay",
     &VectorMesonTensorVectorDecayer::setUpDecayMode, false);

}

void VectorMesonTensorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[1],
					outgoing,true,decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double VectorMesonTensorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin2,PDT::Spin1)));
  // check for photons
  bool photon = outgoing[1]->id()==ParticleID::gamma;
  // // stuff for incoming particle
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel);
    tensors_[ihel] = twave.wave();
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  Energy2 denom[2] = {momenta[0]*part.momentum()-momenta[0].mass()*   part   .mass(),
		      momenta[0]*momenta[1]     -momenta[0].mass()*momenta[1].mass()};
  double fact(coupling_[imode()]/part.mass());
  // calculate the matrix element
  for(unsigned int ih1=0;ih1<5;++ih1) {
    LorentzPolarizationVectorE v1 = tensors_[ih1].preDot (part.momentum());
    LorentzPolarizationVectorE v2 = tensors_[ih1].postDot(momenta[1]);
    complex<Energy2> d0 = v1*momenta[1];
    for(unsigned int ih0=0;ih0<3;++ih0) {
      LorentzPolarizationVector v3 = tensors_[ih1].postDot(vectors_[0][ih0]);
      complex<InvEnergy> d1 = vectors_[0][ih0]*momenta[0]/denom[0];
      complex<Energy> d4 = v2*vectors_[0][ih0];
      for(unsigned int ih2=0;ih2<3;++ih2) {
	complex<InvEnergy> d2 = vectors_[1][ih2]*momenta[0]/denom[1];
	if ( photon && ih2==1)
	  (*ME())(ih0,ih1,ih2)=0.;
	else {
	  (*ME())(ih0,ih1,ih2)=fact*(v3*vectors_[1][ih2] -(v1*vectors_[1][ih2])*d1
				     -d4*d2+d0*d1*d2);
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test = 5.*sqr(coupling_[imode()]/part.mass())/3.;
  // if(photon) test *=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool VectorMesonTensorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
	imode=ix;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
	imode=ix;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
	imode=ix;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
	imode=ix;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=19;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

void VectorMesonTensorVectorDecayer::dataBaseOutput(ofstream & output,
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

string VectorMesonTensorVectorDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
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
