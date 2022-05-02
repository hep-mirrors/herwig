// -*- C++ -*-
//
// TensorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorVectorDecayer class.
//

#include "TensorMesonVectorVectorDecayer.h"
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

void TensorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorVectorDecayer" 
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

int TensorMesonVectorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
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


void TensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV);
}

void TensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigTensorMesonVectorVectorDecayer("Herwig::TensorMesonVectorVectorDecayer", "HwTMDecay.so");

void TensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorVectorDecayer> documentation
    ("The TensorMesonVectorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static Command<TensorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &TensorMesonVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<TensorMesonVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

}
void TensorMesonVectorVectorDecayer::
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
double TensorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  // photons ??
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  for(unsigned int iy=0;iy<2;++iy) {
    vectors_[iy].resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(photon[iy] && ix==1) continue;
      vectors_[iy][ix] = HelicityFunctions::polarizationVector(-momenta[iy],ix,Helicity::outgoing);
    }
  }
  // compute some useful dot products etc
  complex<Energy> p1eps2[3],p2eps1[3];
  Energy2 p1p2(momenta[0]*momenta[1]);
  for(unsigned int ix=0;ix<3;++ix) {
    p1eps2[ix]=vectors_[1][ix]*momenta[0];
    p2eps1[ix]=vectors_[0][ix]*momenta[1];
  }
  // compute the traces and useful dot products
  Complex trace[5];
  complex<Energy2> pboth[5];
  LorentzPolarizationVectorE pleft[2][5],pright[2][5];
  for(unsigned int ix=0;ix<5;++ix) {
    trace[ix]=tensors_[ix].xx() + tensors_[ix].yy() + tensors_[ix].zz() + tensors_[ix].tt();
    pleft[0][ix] =(-momenta[0])*tensors_[ix];
    pleft[1][ix] =(-momenta[1])*tensors_[ix];
    pright[0][ix]=tensors_[ix]*(-momenta[0]);
    pright[1][ix]=tensors_[ix]*(-momenta[1]);
    pboth[ix]=((pleft[0][ix]+pright[0][ix])*momenta[1]);
  }
  // loop to compute the matrix element
  Complex e1e2;
  LorentzTensor<Energy2> te1e2;
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  complex<Energy2> me;
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon[0]&&ix==1) continue;
    for(unsigned iy=0;iy<3;++iy) {
      if(photon[1]&&iy==1) continue;
      e1e2=vectors_[0][ix].dot(vectors_[1][iy]);
      te1e2=complex<Energy2>(p1p2)*
	(LorentzTensor<double>(vectors_[0][ix],vectors_[1][iy])+
	 LorentzTensor<double>(vectors_[1][iy],vectors_[0][ix]));
      for(unsigned int inhel=0;inhel<5;++inhel) {
	me = (tensors_[inhel]*te1e2
	      -p2eps1[ix]*(vectors_[1][iy].dot(pleft[0][inhel]+pright[0][inhel])) 
	      -p1eps2[iy]*(vectors_[0][ix].dot(pleft[1][inhel]+pright[1][inhel]))
	      +pboth[inhel]*e1e2
	      +(p2eps1[ix]*p1eps2[iy]-e1e2*p1p2)*trace[inhel]);
	(*ME())(inhel,ix,iy)=Complex(fact*me);
      }    
    }
  }
  double output = ME()->contract(rho_).real();
  // testing the matrix element 
  // Energy2 m02(sqr(part.mass())),m12(sqr(momenta[0].mass())),m22(sqr(momenta[1].mass()));
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass());
  // Energy2 pcm2(sqr(pcm));
  // double test = 4./15./m02/m02*sqr(coupling_[imode()])*
  //   (3.*m02*(8.*pcm2*pcm2+5.*(m12*m22+pcm2*(m12+m22)))-5.*(m12-m22)*(m12-m22)*pcm2);
  // cout << "testing the matrix element VV " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test <<  " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool TensorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=6;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

void TensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
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

string TensorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
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
