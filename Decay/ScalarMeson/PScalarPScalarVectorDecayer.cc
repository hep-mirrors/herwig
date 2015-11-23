// -*- C++ -*-
//
// PScalarPScalarVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarPScalarVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

PScalarPScalarVectorDecayer::PScalarPScalarVectorDecayer() 
  : _incoming(19), _outgoingP(19), _outgoingV(19), 
    _coupling(19), _maxweight(19) {
  // decay pi' to rho pi
  _incoming[0] =  100111; _outgoingV[0] =  213; _outgoingP[0] = -211; 
  _coupling[0] = 3.57; _maxweight[0] = 4.5; 
  _incoming[1] =  100211; _outgoingV[1] =  213; _outgoingP[1] =  111; 
  _coupling[1] = 3.57; _maxweight[1] = 4.5; 
  _incoming[2] =  100211; _outgoingV[2] =  113; _outgoingP[2] =  211; 
  _coupling[2] = 3.57; _maxweight[2] = 4.5; 
  // K' to K rho
  _incoming[3] =  100311; _outgoingP[3] =  311; _outgoingV[3] =  113; 
  _coupling[3] = 1.; _maxweight[3] = 4.; 
  _incoming[4] =  100321; _outgoingP[4] =  321; _outgoingV[4] =  113; 
  _coupling[4] = 1.; _maxweight[4] = 4.; 
  _incoming[5] =  100311; _outgoingP[5] =  321; _outgoingV[5] = -213; 
  _coupling[5] = 1.41; _maxweight[5] = 4.; 
  _incoming[6] =  100321; _outgoingP[6] =  311; _outgoingV[6] =  213; 
  _coupling[6] = 1.41; _maxweight[6] = 4.; 
  // K' to K* pi
  _incoming[7] =  100311; _outgoingV[7] =  313; _outgoingP[7] =  111; 
  _coupling[7] = 1.55; _maxweight[7] = 2.; 
  _incoming[8] =  100321; _outgoingV[8] =  323; _outgoingP[8] =  111; 
  _coupling[8] = 1.55; _maxweight[8] = 2.; 
  _incoming[9] =  100311; _outgoingV[9] =  323; _outgoingP[9] = -211; 
  _coupling[9] = 2.19; _maxweight[9] = 2.; 
  _incoming[10] =  100321; _outgoingV[10] =  313; _outgoingP[10] =  211; 
  _coupling[10] = 2.19; _maxweight[10] = 2.; 
  // eta (1475) to K* K
  _incoming[11] =  100331; _outgoingV[11] =  323; _outgoingP[11] = -321; 
  _coupling[11] = 2.92; _maxweight[11] = 3.5; 
  _incoming[12] =  100331; _outgoingV[12] =  313; _outgoingP[12] = -311; 
  _coupling[12] = 2.92; _maxweight[12] = 3.5; 
  // eta (1475) to K* K
  _incoming[13] =  9020221; _outgoingV[13] =  323; _outgoingP[13] = -321; 
  _coupling[13] = 0.956; _maxweight[13] = 4.; 
  _incoming[14] =  9020221; _outgoingV[14] =  313; _outgoingP[14] = -311; 
  _coupling[14] = 0.956; _maxweight[14] = 4.; 
  // decay f_0(1370) to a_1 pi
  _incoming[15] =  10221; _outgoingV[15] = 20213; _outgoingP[15] = -211; 
  _coupling[15] = 2.68; _maxweight[15] = 4.5; 
  _incoming[16] =  10221; _outgoingV[16] = 20113; _outgoingP[16] =  111; 
  _coupling[16] = 2.68; _maxweight[16] = 4.5; 
  // decay f_0(1500) to a_1 pi
  _incoming[17] =  9030221; _outgoingV[17] = 20213; _outgoingP[17] = -211; 
  _coupling[17] = 1.147; _maxweight[17] = 3.2; 
  _incoming[18] =  9030221; _outgoingV[18] = 20113; _outgoingP[18] =  111; 
  _coupling[18] = 1.147; _maxweight[18] = 3.2; 
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
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingP[ix]);
    extpart[2] = getParticleData(_outgoingV[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
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

ClassDescription<PScalarPScalarVectorDecayer> 
PScalarPScalarVectorDecayer::initPScalarPScalarVectorDecayer;
// Definition of the static class description member.

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

double PScalarPScalarVectorDecayer::me2( const int,
					 const Particle & inpart,
					 const ParticleVector & decay,
					 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
    VectorWaveFunction::constructSpinInfo(_vectors,decay[1],
					  outgoing,true,false);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors,decay[1],
					     outgoing,false);
  // calculate the matrix element
  Lorentz5Momentum psum(inpart.momentum()+decay[0]->momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(0,0,ix)=_coupling[imode()]/inpart.mass()*(_vectors[ix]*psum);
  }
  // test of the matrix element
//   double me=newME.contract(rhoin).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = 4.*sqr(_coupling[imode()]*pcm/decay[1]->mass());
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return ME()->contract(_rho).real();
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
