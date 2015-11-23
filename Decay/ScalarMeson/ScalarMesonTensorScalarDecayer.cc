// -*- C++ -*-
//
// ScalarMesonTensorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonTensorScalarDecayer class.
//

#include "ScalarMesonTensorScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarMesonTensorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix) 
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

ScalarMesonTensorScalarDecayer::ScalarMesonTensorScalarDecayer() 
  : _incoming(3), _outgoingT(3), _outgoingS(3), _coupling(3), _maxweight(3) {
  // D+ -> f_2 pi
  _incoming[0] =  411; _outgoingT[0] = 225; _outgoingS[0] =  211; 
  _coupling[0] = 8.23E-7/GeV; _maxweight[0] = 5; 
  // chi_c0 -> K*_0 K*_2
  _incoming[1] =  10441; _outgoingT[1] = 325; _outgoingS[1] =  -10321; 
  _coupling[1] = 0.0217/GeV; _maxweight[1] = 5; 
  _incoming[2] =  10441; _outgoingT[2] = 315; _outgoingS[2] =  -10311; 
  _coupling[2] = 0.0217/GeV; _maxweight[2] = 5; 
  // initial size of the arrays
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void ScalarMesonTensorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=_coupling.size();
  if(isize!=_incoming.size()  || isize!=_outgoingT.size()||
     isize!=_outgoingS.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in "
			  << "ScalarMesonTensorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingT[ix]);
    extpart[2] = getParticleData(_outgoingS[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int ScalarMesonTensorScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(id0   ==_incoming[ix]) {
      if((id1   ==_outgoingT[ix]&&id2   ==_outgoingS[ix])||
	 (id2   ==_outgoingT[ix]&&id1   ==_outgoingS[ix])) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==_incoming[ix]&&imode<0) {
      if((id1bar==_outgoingT[ix]&&id2bar==_outgoingS[ix])||
	 (id2bar==_outgoingT[ix]&&id1bar==_outgoingS[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void ScalarMesonTensorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) << _incoming << _outgoingT << _outgoingS << _maxweight;
}

void ScalarMesonTensorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV) >> _incoming >> _outgoingT >> _outgoingS >> _maxweight;
}

ClassDescription<ScalarMesonTensorScalarDecayer> 
ScalarMesonTensorScalarDecayer::initScalarMesonTensorScalarDecayer;
// Definition of the static class description member.

void ScalarMesonTensorScalarDecayer::Init() {

  static ClassDocumentation<ScalarMesonTensorScalarDecayer> documentation
    ("The ScalarMesonTensorScalarDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarMesonTensorScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingT
    ("OutgoingTensor",
     "The PDG code for the outgoing tensor",
     &ScalarMesonTensorScalarDecayer::_outgoingT,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,int> interfaceOutcomingS
    ("OutgoingScalar",
     "The PDG code for the outgoing scalar",
     &ScalarMesonTensorScalarDecayer::_outgoingS,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarMesonTensorScalarDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 100./GeV, false, false, true);

  static ParVector<ScalarMesonTensorScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarMesonTensorScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

double ScalarMesonTensorScalarDecayer::me2(const int,
					   const Particle & inpart,
					   const ParticleVector & decay,
					   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin2,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    TensorWaveFunction::constructSpinInfo(_tensors,decay[0],
					  outgoing,true,false);
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  TensorWaveFunction::
    calculateWaveFunctions(_tensors,decay[0],outgoing,false);
  // calculate the matrix element
  InvEnergy2 fact(_coupling[imode()]/inpart.mass());
  LorentzPolarizationVectorE vtemp;
  for(unsigned int ix=0;ix<5;++ix) {
    vtemp = _tensors[ix]*inpart.momentum(); 
    (*ME())(0,ix,0) = fact * decay[1]->momentum().dot(vtemp);
  }
  // test of the matrix element
//   double me=newME.contract(rhoin).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = 2.*pow<4,1>(pcm)*sqr(_coupling[imode()]*inpart.mass())/
//     3./pow<4,1>(decay[0]->mass());
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return ME()->contract(_rho).real();
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarMesonTensorScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
						   double & coupling) const {
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
  int imode(-1);
  do {
    if(id   ==_incoming[ix]) {
      if(id1==_outgoingT[ix]&&id2==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoingT[ix]&&id1==_outgoingS[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingT[ix]&&id2bar==_outgoingS[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingT[ix]&&id1bar==_outgoingS[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode]*dm.parent()->mass();
  itype = 11;
  return order;
}

// output the setup information for the particle database
void ScalarMesonTensorScalarDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingTensor " << ix << " " 
	     << _outgoingT[ix] << "\n";
      output << "newdef " << name() << ":OutgoingScalar " << ix << " " 
	     << _outgoingS[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingTensor " << ix << " " 
	     << _outgoingT[ix] << "\n";
      output << "insert " << name() << ":OutgoingScalar " << ix << " " 
	     << _outgoingS[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
