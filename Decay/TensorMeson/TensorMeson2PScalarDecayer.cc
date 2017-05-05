// -*- C++ -*-
//
// TensorMeson2PScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMeson2PScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

TensorMeson2PScalarDecayer::TensorMeson2PScalarDecayer() 
  : _incoming(48), _outgoing1(48), _outgoing2(48), 
    _coupling(48), _maxweight(48) {
  // a_2 -> eta pi
  _incoming[0] = 115; _outgoing1[0] =  221; _outgoing2[0] = 111; 
  _coupling[0] = 10.90/GeV; _maxweight[0] = 1.7; 
  _incoming[1] = 215; _outgoing1[1] =  221; _outgoing2[1] = 211; 
  _coupling[1] = 10.90/GeV; _maxweight[1] = 1.7; 
  // a_2 -> eta' pi
  _incoming[2] = 115; _outgoing1[2] =  331; _outgoing2[2] = 111; 
  _coupling[2] = 9.92/GeV; _maxweight[2] = 1.9; 
  _incoming[3] = 215; _outgoing1[3] =  331; _outgoing2[3] = 211; 
  _coupling[3] = 9.92/GeV; _maxweight[3] = 1.9; 
  // a_2 -> K K
  _incoming[4] = 115; _outgoing1[4] =  311; _outgoing2[4] = -311; 
  _coupling[4] = 7.36/GeV; _maxweight[4] = 1.7; 
  _incoming[5] = 115; _outgoing1[5] =  321; _outgoing2[5] = -321; 
  _coupling[5] = 7.36/GeV; _maxweight[5] = 1.7; 
  _incoming[6] = 215; _outgoing1[6] =  321; _outgoing2[6] = -311; 
  _coupling[6] = 10.41/GeV; _maxweight[6] = 1.7; 
  // f_2 -> pi pi
  _incoming[7] = 225; _outgoing1[7] =  211; _outgoing2[7] = -211; 
  _coupling[7] = 18.73/GeV; _maxweight[7] = 1.7; 
  _incoming[8] = 225; _outgoing1[8] =  111; _outgoing2[8] = 111; 
  _coupling[8] = 13.24/GeV; _maxweight[8] = 1.7; 
  // f_2 -> eta eta
  _incoming[9] = 225; _outgoing1[9] =  221; _outgoing2[9] = 221; 
  _coupling[9] = 8.362/GeV; _maxweight[9] = 1.8; 
  // f_2 -> KK
  _incoming[10] = 225; _outgoing1[10] =  321; _outgoing2[10] = -321; 
  _coupling[10] = 11.03/GeV; _maxweight[10] = 2.; 
  _incoming[11] = 225; _outgoing1[11] =  311; _outgoing2[11] = -311; 
  _coupling[11] = 11.38/GeV; _maxweight[11] = 2.; 
  // f'_2 -> KK
  _incoming[12] = 335; _outgoing1[12] =  321; _outgoing2[12] = -321; 
  _coupling[12] = 14.65/GeV; _maxweight[12] = 1.65; 
  _incoming[13] = 335; _outgoing1[13] =  311; _outgoing2[13] = -311; 
  _coupling[13] = 14.65/GeV; _maxweight[13] = 1.65; 
  // f'_2 -> eta eta
  _incoming[14] = 335; _outgoing1[14] =  221; _outgoing2[14] = 221; 
  _coupling[14] = 9.15/GeV; _maxweight[14] = 1.7; 
  // f'_2 -> pi pi 
  _incoming[15] = 335; _outgoing1[15] =  211; _outgoing2[15] = -211; 
  _coupling[15] = 0.860/GeV; _maxweight[15] = 1.7; 
  _incoming[16] = 335; _outgoing1[16] =  111; _outgoing2[16] = 111; 
  _coupling[16] = 0.608/GeV; _maxweight[16] = 1.7; 
  // K_2 -> K eta
  _incoming[17] =  325; _outgoing1[17] =  321; _outgoing2[17] = 221; 
  _coupling[17] = 1.52/GeV; _maxweight[17] = 1.8; 
  _incoming[18] =  315; _outgoing1[18] =  311; _outgoing2[18] = 221; 
  _coupling[18] = 1.52/GeV; _maxweight[18] = 1.8; 
  // K_2 -> K pi
  _incoming[19] =  325; _outgoing1[19] =  321; _outgoing2[19] =  111; 
  _coupling[19] = 8.30/GeV;_maxweight[19] = 1.65; 
  _incoming[20] =  325; _outgoing1[20] =  311; _outgoing2[20] =  211; 
  _coupling[20] = 11.74/GeV; _maxweight[20] = 1.65; 
  _incoming[21] =  315; _outgoing1[21] =  311; _outgoing2[21] =  111; 
  _coupling[21] = 8.68/GeV; _maxweight[21] = 1.65; 
  _incoming[22] =  315; _outgoing1[22] =  321; _outgoing2[22] = -211; 
  _coupling[22] = 12.28/GeV; _maxweight[22] = 1.65; 
  // B_2 -> B pi
  _incoming[23] =  525; _outgoing1[23] =  521; _outgoing2[23] = 111; 
  _coupling[23] = 27.23/GeV; _maxweight[23] = 1.65; 
  _incoming[24] =  525; _outgoing1[24] =  511; _outgoing2[24] = 211; 
  _coupling[24] = 38.52/GeV; _maxweight[24] = 1.65; 
  _incoming[25] =  515; _outgoing1[25] =  511; _outgoing2[25] = 111; 
  _coupling[25] = 27.16/GeV; _maxweight[25] = 1.65; 
  _incoming[26] =  515; _outgoing1[26] =  521; _outgoing2[26] = -211; 
  _coupling[26] = 38.62/GeV; _maxweight[26] = 1.65; 
  // D_2 -> D pi
  _incoming[27] =  425; _outgoing1[27] =  421; _outgoing2[27] = 111; 
  _coupling[27] = 18.07/GeV; _maxweight[27] = 1.7; 
  _incoming[28] =  425; _outgoing1[28] =  411; _outgoing2[28] = -211; 
  _coupling[28] = 25.56/GeV; _maxweight[28] = 1.7; 
  _incoming[29] =  415; _outgoing1[29] =  411; _outgoing2[29] = 111; 
  _coupling[29] = 14.91/GeV; _maxweight[29] = 1.7; 
  _incoming[30] =  415; _outgoing1[30] =  421; _outgoing2[30] = 211; 
  _coupling[30] = 21.09/GeV; _maxweight[30] = 1.7; 
  // D_s2
  _incoming[31] =  435; _outgoing1[31] =  421; _outgoing2[31] =  321; 
  _coupling[31] = 23.39/GeV; _maxweight[31] = 1.7; 
  _incoming[32] =  435; _outgoing1[32] =  411; _outgoing2[32] =  311; 
  _coupling[32] = 23.39/GeV; _maxweight[32] = 1.7; 
  // B_s2
  _incoming[33] =  535; _outgoing1[33] =  521; _outgoing2[33] = -321; 
  _coupling[33] = 45.43/GeV; _maxweight[33] = 1.7; 
  _incoming[34] =  535; _outgoing1[34] =  511; _outgoing2[34] = -311; 
  _coupling[34] = 48.84/GeV; _maxweight[34] = 1.7; 
  // chi_c2 to pi pi
  _incoming[35] =  445; _outgoing1[35] =  211; _outgoing2[35] = -211; 
  _coupling[35] = 0.0200/GeV; _maxweight[35] = 1.7; 
  _incoming[36] =  445; _outgoing1[36] =  111; _outgoing2[36] =  111; 
  _coupling[36] = 0.0141/GeV; _maxweight[36] = 1.7; 
  // chi_c2 to K K
  _incoming[37] =  445; _outgoing1[37] =  321; _outgoing2[37] = -321; 
  _coupling[37] = 0.056/GeV; _maxweight[37] = 1.7; 
  _incoming[38] =  445; _outgoing1[38] =  311; _outgoing2[38] = -311; 
  _coupling[38] = 0.056/GeV; _maxweight[38] = 1.7; 
  // f_2 to sigma sigma
  _incoming[39] = 225; _outgoing1[39] = 9000221; _outgoing2[39] = 9000221; 
  _coupling[39] = 104.1/GeV; _maxweight[39] = 60.; 
  // pi_2 to sigma pi
  _incoming[40] = 10115; _outgoing1[40] = 9000221; _outgoing2[40] = 111; 
  _coupling[40] = 15.3/GeV; _maxweight[40] = 10.; 
  _incoming[41] = 10215; _outgoing1[41] = 9000221; _outgoing2[41] = 211; 
  _coupling[41] = 15.3/GeV; _maxweight[41] = 10.; 
  // eta_2 to a_0 pi
  _incoming[42] = 10225; _outgoing1[42] = 9000111; _outgoing2[42] = 111; 
  _coupling[42] = 11.3/GeV; _maxweight[42] = 10.; 
  _incoming[43] = 10225; _outgoing1[43] = 9000211; _outgoing2[43] = -211; 
  _coupling[43] = 11.3/GeV; _maxweight[43] = 10.; 
  // eta'_2 to a_0 pi
  _incoming[44] = 10335; _outgoing1[44] = 9000111; _outgoing2[44] = 111; 
  _coupling[44] = 4.43/GeV; _maxweight[44] = 10.; 
  _incoming[45] = 10335; _outgoing1[45] = 9000211; _outgoing2[45] = -211; 
  _coupling[45] = 4.43/GeV; _maxweight[45] = 10.; 
  // chi_c2(2P) to D D
  _incoming[46] =  100445; _outgoing1[46] =  411; _outgoing2[46] = -411; 
  _coupling[46] = 22.3/GeV; _maxweight[46] = 1.7; 
  _incoming[47] =  100445; _outgoing1[47] =  421; _outgoing2[47] = -421; 
  _coupling[47] = 22.3/GeV; _maxweight[47] = 1.7; 
  // initial size of the vectors
  _initsize=_incoming.size();
  // intermediates
  generateIntermediates(false);
}

void TensorMeson2PScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters TensorMeson2PScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoing1[ix]);
    extpart[2]=getParticleData(_outgoing2[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int TensorMeson2PScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
      if((id1   ==_outgoing1[ix]&&id2   ==_outgoing2[ix])||
	 (id2   ==_outgoing1[ix]&&id1   ==_outgoing2[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
      if((id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix])||
	 (id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}

void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << _incoming << _outgoing1 << _outgoing2 << _maxweight 
     << ounit(_coupling,1/GeV);
}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight 
     >> iunit(_coupling,1/GeV);
}

ClassDescription<TensorMeson2PScalarDecayer> 
TensorMeson2PScalarDecayer::initTensorMeson2PScalarDecayer;
// Definition of the static class description member.

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The TensorMeson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMeson2PScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMeson2PScalarDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMeson2PScalarDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 1000./GeV, false, false, true);

  static ParVector<TensorMeson2PScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMeson2PScalarDecayer::_maxweight,
     0, 0, 0, 0., 100000., false, false, true);
}

// matrix elememt for the process
double TensorMeson2PScalarDecayer::me2(const int, const Particle & inpart,
				       const ParticleVector & decay,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
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
    for(unsigned int ix=0;ix<decay.size();++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return 0.;
  }
  // calculate the matrix element
  for(unsigned int ix=0;ix<5;++ix) {
    (*ME())(ix,0,0) = _coupling[imode()]/inpart.mass()*
      ((_tensors[ix]*decay[1]->momentum())*decay[0]->momentum());
  }
//   // test of the answer
//   double me = newME.contract(_rho).real();
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   double test = Energy4(pow<4,1>(2*pcm))*sqr( _coupling[imode()]/inpart.mass())/120.;
//   cout << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return ME()->contract(_rho).real();
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
    if(id   ==_incoming[ix]) {
      if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=7;
  return order;
}

void TensorMeson2PScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "newdef " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":FirstOutgoing " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << name() << ":SecondOutgoing " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix]*GeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
