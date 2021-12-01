// -*- C++ -*-
//
// VectorMeson2FermionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2FermionDecayer class.
//

#include "VectorMeson2FermionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2FermionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2FermionDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;
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

VectorMeson2FermionDecayer::VectorMeson2FermionDecayer() {
  // don't include intermediates
  generateIntermediates(false);
}

int VectorMeson2FermionDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(incoming_[ix]==id   ) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(incoming_[ix]==idbar) {
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

void VectorMeson2FermionDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2FermionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMeson2FermionDecayer,DecayIntegrator>
describeHerwigVectorMeson2FermionDecayer("Herwig::VectorMeson2FermionDecayer", "HwVMDecay.so");

void VectorMeson2FermionDecayer::Init() {

  static ClassDocumentation<VectorMeson2FermionDecayer> documentation
    ("The VectorMeson2FermionDecayer class is designed for the decay "
     "of vector mesons to fermions. It is mainly used for the decay of vector mesons "
     "to electrons and muons.");
  
  static Command<VectorMeson2FermionDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, fermion, antifermion), coupling and max weight for a decay",
     &VectorMeson2FermionDecayer::setUpDecayMode, false);

  static Deleted<VectorMeson2FermionDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceOutcoming1
    ("OutgoingFermion","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceOutcoming2
    ("OutgoingAntiFermion","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMeson2FermionDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double VectorMeson2FermionDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // prefactor
  InvEnergy pre(coupling_[imode()]/part.mass());
  // now compute the currents
  LorentzPolarizationVector temp;
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      temp = pre*wave_[ix].vectorCurrent(wavebar_[iy]);
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 4.*sqr(coupling_[imode()])/3.*
  //   (1.+2.*sqr(momenta[0].mass()/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMeson2FermionDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
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
  mecode=2;
  return order;
}

// output the setup information for the particle database
void VectorMeson2FermionDecayer::dataBaseOutput(ofstream & output,
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
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMeson2FermionDecayer::setUpDecayMode(string arg) {
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
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
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
