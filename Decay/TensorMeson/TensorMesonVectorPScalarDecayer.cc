// -*- C++ -*-
//
// TensorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorPScalarDecayer class.
//

#include "TensorMesonVectorPScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

TensorMesonVectorPScalarDecayer::TensorMesonVectorPScalarDecayer() 
  :  _incoming(31), _outgoingV(31), _outgoingP(31), 
     _coupling(31), _maxweight(31) {
  // a_2 -> rho pi
  _incoming[0] =  115; _outgoingV[0] =  213; _outgoingP[0] = -211; 
  _coupling[0] = 21.1/GeV2; _maxweight[0] = 10.; 
  _incoming[1] =  215; _outgoingV[1] =  113; _outgoingP[1] =  211; 
  _coupling[1] = 21.1/GeV2; _maxweight[1] = 9.; 
  _incoming[2] =  215; _outgoingV[2] =  213; _outgoingP[2] =  111; 
  _coupling[2] = 21.1/GeV2; _maxweight[2] = 9.; 
  // a_2+/- -> gamma pi+/-
  _incoming[3] =  215; _outgoingV[3] =  22; _outgoingP[3] =  211; 
  _coupling[3] = 0.551/GeV2; _maxweight[3] = 2.; 
  // k_2 -> K_2 omega
  _incoming[4] =  315; _outgoingV[4] = 223; _outgoingP[4] =  311; 
  _coupling[4] = 11.66/GeV2; _maxweight[4] = 17.; 
  _incoming[5] =  325; _outgoingV[5] = 223; _outgoingP[5] =  321; 
  _coupling[5] = 11.66/GeV2; _maxweight[5] = 20.5; 
  // k_2+/- -> K+/- gamma
  _incoming[6] =  325; _outgoingV[6] = 22; _outgoingP[6] =  321; 
  _coupling[6] = 0.553/GeV2; _maxweight[6] = 2.2; 
  // B_c2 -> B_c gamma
  _incoming[7] =  545; _outgoingV[7] = 22; _outgoingP[7] =  541; 
  _coupling[7] = 0.651/GeV2; _maxweight[7] = 2.; 
  // K_2 -> K rho
  _incoming[8] =  325; _outgoingV[8] =  113; _outgoingP[8] =  321; 
  _coupling[8] = 10.14/GeV2; _maxweight[8] = 9.; 
  _incoming[9] =  325; _outgoingV[9] =  213; _outgoingP[9] =  311; 
  _coupling[9] = 14.33/GeV2; _maxweight[9] = 9.; 
  _incoming[10] =  315; _outgoingV[10] =  113; _outgoingP[10] =  311; 
  _coupling[10] = 10.14/GeV2; _maxweight[10] = 9.; 
  _incoming[11] =  315; _outgoingV[11] = -213; _outgoingP[11] =  321; 
  _coupling[11] = 14.33/GeV2; _maxweight[11] = 9.; 
  // K_2 -> K* pi 
  _incoming[12] =  325; _outgoingV[12] =  323; _outgoingP[12] =  111; 
  _coupling[12] = 9.733/GeV2; _maxweight[12] = 13; 
  _incoming[13] =  325; _outgoingV[13] =  313; _outgoingP[13] =  211; 
  _coupling[13] = 13.77/GeV2; _maxweight[13] = 11; 
  _incoming[14] =  315; _outgoingV[14] =  313; _outgoingP[14] =  111; 
  _coupling[14] = 9.733/GeV2; _maxweight[14] = 8.; 
  _incoming[15] =  315; _outgoingV[15] =  323; _outgoingP[15] = -211; 
  _coupling[15] = 13.77/GeV2; _maxweight[15] = 8.; 
  // D_2 -> D* pi 
  _incoming[16] =  425; _outgoingV[16] =  423; _outgoingP[16] =  111; 
  _coupling[16] = 8.035/GeV2; _maxweight[16] = 2.2; 
  _incoming[17] =  425; _outgoingV[17] =  413; _outgoingP[17] = -211; 
  _coupling[17] = 11.670/GeV2; _maxweight[17] = 2.4; 
  _incoming[18] =  415; _outgoingV[18] =  413; _outgoingP[18] =  111; 
  _coupling[18] = 6.801/GeV2; _maxweight[18] = 2.4; 
  _incoming[19] =  415; _outgoingV[19] =  423; _outgoingP[19] =  211; 
  _coupling[19] = 9.527/GeV2; _maxweight[19] = 2.; 
  // D_s2 -> D* K
  _incoming[20] =  435; _outgoingV[20] =  423; _outgoingP[20] =  321; 
  _coupling[20] = 13.10/GeV2; _maxweight[20] = 2.2; 
  _incoming[21] =  435; _outgoingV[21] =  413; _outgoingP[21] =  311; 
  _coupling[21] = 13.10/GeV2; _maxweight[21] = 2.5; 
  // B_2 -> B* pi 
  _incoming[22] =  525; _outgoingV[22] =  523; _outgoingP[22] =  111; 
  _coupling[22] = 4.99/GeV2; _maxweight[22] = 2.1; 
  _incoming[23] =  525; _outgoingV[23] =  513; _outgoingP[23] =  211; 
  _coupling[23] = 7.059/GeV2; _maxweight[23] = 2.1; 
  _incoming[24] =  515; _outgoingV[24] =  513; _outgoingP[24] =  111; 
  _coupling[24] = 4.99/GeV2; _maxweight[24] = 2.1; 
  _incoming[25] =  515; _outgoingV[25] =  523; _outgoingP[25] = -211; 
  _coupling[25] = 7.059/GeV2; _maxweight[25] = 2.1; 
  // D_s2
  _incoming[26] =  435; _outgoingV[26] =  423; _outgoingP[26] =  321; 
  _coupling[26] = 13.09/GeV2; _maxweight[26] = 2.2; 
  _incoming[27] =  435; _outgoingV[27] =  413; _outgoingP[27] =  311; 
  _coupling[27] = 13.09/GeV2; _maxweight[27] = 2.5; 
  // B_s2
  _incoming[28] =  535; _outgoingV[28] =  523; _outgoingP[28] = -321; 
  _coupling[28] = 7.29/GeV2; _maxweight[28] = 2.4; 
  _incoming[29] =  535; _outgoingV[29] =  513; _outgoingP[29] = -311; 
  _coupling[29] = 9.43/GeV2; _maxweight[29] = 2.1; 
  // upsilon_2(1d) to chi_b gamma
  _incoming[30] = 20555; _outgoingV[30] =   22; _outgoingP[30] = 10551; 
  _coupling[30] = 1.11/GeV2; _maxweight[30] = 2.4; 
  // initial size of the arrays
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void TensorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingV.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorPScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoingV[ix]);
    extpart[2] = getParticleData(_outgoingP[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int TensorMesonVectorPScalarDecayer::modeNumber(bool & cc, tcPDPtr parent, 
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

void TensorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoingV << _outgoingP << _maxweight << ounit(_coupling,1/GeV2);
}

void TensorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoingV >> _outgoingP >> _maxweight >> iunit(_coupling,1/GeV2);
}

ClassDescription<TensorMesonVectorPScalarDecayer> 
TensorMesonVectorPScalarDecayer::initTensorMesonVectorPScalarDecayer;
// Definition of the static class description member.

void TensorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorPScalarDecayer> documentation
    ("The TensorMesonVectorPScalarDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingV
    ("OutgoingVector",
     "The PDG code for the outgoing spin-1particle",
     &TensorMesonVectorPScalarDecayer::_outgoingV,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,int> interfaceOutcomingP
    ("OutgoingScalar",
     "The PDG code for the outgoing pseudoscalar meson",
     &TensorMesonVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,InvEnergy2> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorPScalarDecayer::_coupling,
     1/GeV2, 0, ZERO, ZERO, 100./GeV2, false, false, true);

  static ParVector<TensorMesonVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 1000., false, false, true);

}

// matrix elememt for the process
double TensorMesonVectorPScalarDecayer::me2(const int,const Particle & inpart,
					    const ParticleVector & decay,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin0)));
  // check for photons
  bool photon(_outgoingV[imode()]==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    _rho = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(_tensors,_rho,const_ptr_cast<tPPtr>(&inpart),
			     incoming,false);
  }
  if(meopt==Terminate) {
    TensorWaveFunction::constructSpinInfo(_tensors,const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    // set up the spin information for the decay products
    VectorWaveFunction::constructSpinInfo(_vectors,decay[0],outgoing,true,photon);
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors,decay[0],outgoing,photon);
  InvEnergy3 fact(_coupling[imode()]/inpart.mass());
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    for(unsigned int vhel=0;vhel<3;++vhel){
      if(vhel==1&&photon) (*ME())(inhel,vhel,0)=0.;
      else {
	LorentzVector<complex<InvEnergy> > vtemp=
	  fact*epsilon(decay[0]->momentum(),_vectors[vhel],decay[1]->momentum());
	(*ME())(inhel,vhel,0)= Complex((decay[1]->momentum()*_tensors[inhel]).dot(vtemp));
      }
    }
  }
//   // test of the answer
//   double me = ME()->contract(_rho).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = Energy4(pow<4,1>(2*pcm))*sqr( _coupling[imode()])/80.;
//   cout << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return ME()->contract(_rho).real();
}

bool TensorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling=_coupling[imode]*sqr(dm.parent()->mass());
  mecode=8;
  return order;
}

void TensorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingVector " << ix << " " 
	     << _outgoingV[ix] << "\n";
      output << "newdef " << name() << ":OutgoingScalar " << ix << " " 
	     << _outgoingP[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV2 << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingVector " << ix << " " 
	     << _outgoingV[ix] << "\n";
      output << "insert " << name() << ":OutgoingScalar " << ix << " " 
	     << _outgoingP[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV2 << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
