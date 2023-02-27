// -*- C++ -*-
//
// Spin3MesonTensorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3MesonTensorVectorDecayer class.
//

#include "Spin3MesonTensorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3MesonTensorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3MesonTensorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3MesonTensorVectorDecayer" 
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

int Spin3MesonTensorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void Spin3MesonTensorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void Spin3MesonTensorVectorDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3MesonTensorVectorDecayer,DecayIntegrator>
describeHerwigSpin3MesonTensorVectorDecayer("Herwig::Spin3MesonTensorVectorDecayer", "HwTMDecay.so");

void Spin3MesonTensorVectorDecayer::Init() {

  static ClassDocumentation<Spin3MesonTensorVectorDecayer> documentation
    ("The Spin3MesonTensorVectorDecayer class is designed for the decay"
     " of a tensor meson to a tensor and pseudoscalar mesons.");

    static Command<Spin3MesonTensorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(GeV) and max weight for a decay",
     &Spin3MesonTensorVectorDecayer::setUpDecayMode, false);

}

string Spin3MesonTensorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
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
  out.first = stoi(stype);
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
  Energy g = stof(stype)*GeV;
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

void Spin3MesonTensorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(rank3_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],outgoing,true,
					decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double Spin3MesonTensorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin2,PDT::Spin1)));
  // check for photons
  bool photon(outgoing[1]->id()==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(rank3_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  tensors_.resize(5);
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel,tensor_phase);
    tensors_[ihel] = twave.wave();
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Energy2 denom[2] = {part.momentum()*momenta[0]-part.mass()*momenta[0].mass(),
		      part.momentum()*momenta[1]-part.mass()*momenta[1].mass()};
  double fact(coupling_[imode()]/part.mass());
  for(unsigned int ih0=0;ih0<7;++ih0) {
    for(unsigned int ih2=0;ih2<3;++ih2) {
      if(ih2==1 && photon) {
	for(unsigned int ih1=0;ih1<5;++ih1)  (*ME())(ih0,ih1,ih2)=0.;
      }
      else {
	LorentzPolarizationVector v0  = vectors_[ih2] - momenta[1]*momenta[0].dot(vectors_[ih2])/denom[1];
	LorentzTensor<double>     t0  = rank3_[ih0].dot(v0,0);
	LorentzPolarizationVectorE v1 = t0.postDot(momenta[0]);
	Complex d1 = v1*momenta[0]/denom[0];
	for(unsigned int ih1=0;ih1<5;++ih1) {
	  LorentzPolarizationVectorE v2 = tensors_[ih1].postDot(momenta[1]);
	  Complex d2 = v2*momenta[1]/denom[0];
	  (*ME())(ih0,ih1,ih2) = fact*(t0*tensors_[ih1] - 2.*v1.dot(v2)/denom[0]+d1*d2);
	}
      }
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // double test = sqr(fact);
  // if(photon) test *=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool Spin3MesonTensorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=21;
  return order;
}

void Spin3MesonTensorVectorDecayer::dataBaseOutput(ofstream & output,
						   bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
