// -*- C++ -*-
//
// VectorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorVectorDecayer class.
//

#include "VectorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxweight_[ix]=mode(ix)->maxWeight();
  }
}

void VectorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorVectorDecayer" << Exception::runerror;
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

VectorMesonVectorVectorDecayer::VectorMesonVectorVectorDecayer() {
  // intermediates
  generateIntermediates(false);
}

int VectorMesonVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void VectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxweight_ << coupling_;
}

void VectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxweight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigVectorMesonVectorVectorDecayer("Herwig::VectorMesonVectorVectorDecayer", "HwVMDecay.so");

void VectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorVectorDecayer> documentation
    ("The VectorMesonVectorVectorDecayer class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static Command<VectorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling and max weight for a decay",
     &VectorMesonVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceOutgoing1
    ("Outgoing1","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceOutgoing2
    ("Outgoing2","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

}

string VectorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
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
  coupling_.push_back(g);
  maxweight_.push_back(wgt);
  // success
  return "";
}

void VectorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix+1],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double VectorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix) 
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix+1].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix+1][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],
								   ihel,Helicity::outgoing);
    }
  }
  // work out the dot products we need for the matrix element
  Energy2 p1p2((momenta[0])*(momenta[1]));
  complex<Energy> p1eps2[3],p2eps1[3];
  for(unsigned int ix=0;ix<3;++ix) {
    p1eps2[ix]=vectors_[2][ix]*(momenta[0]);
    p2eps1[ix]=vectors_[1][ix]*(momenta[1]);
  }
  // compute the matrix element
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Energy2 m12(momenta[0].mass()*momenta[0].mass()),m22(momenta[1].mass()*momenta[1].mass());
  InvEnergy3 fact(2.*coupling_[imode()]/(part.mass()*part.mass()*part.mass()));
  LorentzPolarizationVector vtemp;
  for(unsigned int ipol1=0;ipol1<3;++ipol1) {
    for(unsigned int ipol2=0;ipol2<3;++ipol2) {
      Complex eps1eps2=vectors_[1][ipol1].dot(vectors_[2][ipol2]);
      vtemp=fact*(p1eps2[ipol2]*p2eps1[ipol1]*pdiff
		  +p1eps2[ipol2]*m22*vectors_[1][ipol1]
		  -p2eps1[ipol1]*m12*vectors_[2][ipol2]
		  +eps1eps2*(-p1p2*pdiff+m12*momenta[1]
			     -m22*momenta[0]));
      for(unsigned int inpol=0;inpol<3;++inpol) 
	(*ME())(inpol,ipol1,ipol2)=vectors_[0][inpol].dot(vtemp);
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element;
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = 8./3.*sqr(coupling_[imode()]*pcm/part.mass())*
  //   (1.+sqr(momenta[0].mass()/part.mass())+sqr(momenta[1].mass()/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling = coupling_[imode]; 
  mecode = 5;
  return order; 
}

// output the setup information for the particle database
void VectorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header){output << "update decayers set parameters=\"";}
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
