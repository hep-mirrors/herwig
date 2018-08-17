// -*- C++ -*-
//
// PScalarPScalarVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarPScalarVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

PScalarPScalarVectorDecayer::PScalarPScalarVectorDecayer() {
  _incoming  = { 100111, 100211, 100211, 100311, 100321, 100311, 100321,
		 100311, 100321, 100311, 100321, 100331, 100331,9020221,9020221,
	 	  10221,  10221,9030221,9030221};
  _outgoingP = {   -211,    111,    211,    311,    321,    321,    311,
		    111,    111,   -211,    211,   -321,   -311,   -321,   -311,
		   -211,    111,   -211,    111};
  _outgoingV = {    213,    213,    113,    113,    113,   -213,    213,
		    313,    323,    323,    313,    323,    313,    323,    313,
		  20213,  20113,  20213,  20113};
  _coupling  = {   3.57,   3.57,   3.57,     1.,     1.,   1.41,   1.41,
		   1.55,   1.55,   2.19,   2.19,   2.92,   2.92,  0.956,  0.956,
		   2.68,   2.68,  1.147,  1.147};
  _maxweight = {    4.5,    4.5,    4.5,     4.,     4.,     4.,     4.,
		     2.,     2.,     2.,     2.,    3.5,    3.5,     4.,     4.,
		    4.5,    4.5,    3.2,    3.2};
  // initial size of the arrays
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void PScalarPScalarVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingP.size()||
     isize!=_outgoingV.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in PScalarPScalarVectorDecayer" 
  			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    tPDPtr    in  =  getParticleData( _incoming[ix]);
    tPDVector out = {getParticleData(_outgoingP[ix]),
		     getParticleData(_outgoingV[ix])};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,_maxweight[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PScalarPScalarVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(id   ==_incoming[ix]) {
      if((id1   ==_outgoingP[ix]&&id2   ==_outgoingV[ix])||
	 (id2   ==_outgoingP[ix]&&id1   ==_outgoingV[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
      if((id1bar==_outgoingP[ix]&&id2bar==_outgoingV[ix])||
	 (id2bar==_outgoingP[ix]&&id1bar==_outgoingV[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void PScalarPScalarVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _coupling << _incoming << _outgoingP << _outgoingV << _maxweight;
}

void PScalarPScalarVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _coupling >> _incoming >> _outgoingP >> _outgoingV >> _maxweight;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarPScalarVectorDecayer,DecayIntegrator>
describeHerwigPScalarPScalarVectorDecayer("Herwig::PScalarPScalarVectorDecayer", "HwSMDecay.so");

void PScalarPScalarVectorDecayer::Init() {

  static ClassDocumentation<PScalarPScalarVectorDecayer> documentation
    ("The PScalarPScalarVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarPScalarVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &PScalarPScalarVectorDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,int> interfaceOutgoingVector
    ("OutgoingVector",
     "The PDG code for the outgoing vector meson",
     &PScalarPScalarVectorDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarPScalarVectorDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<PScalarPScalarVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarPScalarVectorDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

void PScalarPScalarVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(_vectors,decay[1],
					outgoing,true,false);
}

double PScalarPScalarVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  if(meopt==Initialize) {
     ScalarWaveFunction::
       calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  _vectors.resize(3);
  bool massless = outgoing[1]->id()==ParticleID::gamma;
  for(unsigned int ix=0;ix<3;++ix) {
    if(massless && ix==1) continue;
    _vectors[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Lorentz5Momentum psum(part.momentum()+momenta[0]);
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(0,0,ix)=_coupling[imode()]/part.mass()*(_vectors[ix]*psum);
  }
  double me=ME()->contract(_rho).real();
  // test of the matrix element
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 4.*sqr(_coupling[imode()]*pcm/momenta[1].mass());
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarPScalarVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
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
  unsigned int ix(0); bool order(true);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoingP[ix]&&id2==_outgoingV[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoingP[ix]&&id1==_outgoingV[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingP[ix]&&id2bar==_outgoingV[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingP[ix]&&id1bar==_outgoingV[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode];
  mecode=10;
  return order;
}

// output the setup information for the particle database
void PScalarPScalarVectorDecayer::dataBaseOutput(ofstream & output,
						 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingPScalar " << ix << " " 
	     << _outgoingP[ix] << "\n";
      output << "newdef " << name() << ":OutgoingVector " << ix << " " 
	     << _outgoingV[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingPScalar " << ix << " " 
	     << _outgoingP[ix] << "\n";
      output << "insert " << name() << ":OutgoingVector " << ix << " " 
	     << _outgoingV[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
