// -*- C++ -*-
//
// TensorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorVectorDecayer class.
//

#include "TensorMesonVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
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

TensorMesonVectorVectorDecayer::TensorMesonVectorVectorDecayer() 
  : _incoming(22), _outgoing1(22), _outgoing2(22), 
    _coupling(22), _maxweight(22) {
  // a_2 -> gamma gamma
  _incoming[0] = 115; _outgoing1[0] =  22; _outgoing2[0] = 22; 
  _coupling[0] = 0.00727/GeV; _maxweight[0] = 1.7; 
  // f_2 -> gamma gamma
  _incoming[1] = 225; _outgoing1[1] =  22; _outgoing2[1] = 22; 
  _coupling[1] = 0.01253/GeV; _maxweight[1] = 1.7; 
  // f'_2 -> gamma gamma
  _incoming[2] = 335; _outgoing1[2] =  22; _outgoing2[2] = 22; 
  _coupling[2] = 0.00161/GeV; _maxweight[2] = 1.7; 
  // chi_b(2P) decays
  _incoming[3] = 100555; _outgoing1[3] = 553; _outgoing2[3] = 223; 
  _coupling[3] = 0.0118/GeV; _maxweight[3] = 1.8; 
  _incoming[4] = 100555; _outgoing1[4] = 553; _outgoing2[4] = 22; 
  _coupling[4] = 0.0172/GeV; _maxweight[4] = 1.7; 
  _incoming[5] = 100555; _outgoing1[5] = 100553; _outgoing2[5] = 22; 
  _coupling[5] = 0.145/GeV; _maxweight[5] = 1.7; 
  _incoming[6] = 100555; _outgoing1[6] = 333; _outgoing2[6] = 333; 
  _coupling[6] = 0.00483/GeV; _maxweight[6] = 18.0; 
  // chi_c decays
  _incoming[7] = 445; _outgoing1[7] = 443; _outgoing2[7] = 22; 
  _coupling[7] = 0.243/GeV; _maxweight[7] = 1.7; 
  _incoming[8] = 445; _outgoing1[8] = 323; _outgoing2[8] = -323; 
  _coupling[8] = 0.00560/GeV; _maxweight[8] = 15.; 
  _incoming[9] = 445; _outgoing1[9] = 313; _outgoing2[9] = -313; 
  _coupling[9] = 0.00560/GeV; _maxweight[9] = 20.; 
  _incoming[10] = 445; _outgoing1[10] = 333; _outgoing2[10] = 333; 
  _coupling[10] = 0.00418/GeV; _maxweight[10] = 10.; 
  _incoming[11] = 445; _outgoing1[11] = 22; _outgoing2[11] = 22; 
  _coupling[11] = 0.00122/GeV; _maxweight[11] = 1.7; 
  // chi_b(1P) decays
  _incoming[12] = 555; _outgoing1[12] = 553; _outgoing2[12] = 22; 
  _coupling[12] = 0.0683/GeV; _maxweight[12] = 1.8; 
  // a_2 omega rho 
  _incoming[13] =  115; _outgoing1[13] =  223; _outgoing2[13] =  113; 
  _coupling[13] = 23.1/GeV; _maxweight[13] = 15.; 
  _incoming[14] =  215; _outgoing1[14] =  223; _outgoing2[14] =  213; 
  _coupling[14] = 23.1/GeV; _maxweight[14] = 21.; 
  // f_2 rho rho
  _incoming[15] =  225; _outgoing1[15] =  113; _outgoing2[15] =  113; 
  _coupling[15] = 11.7/GeV; _maxweight[15] = 26.; 
  _incoming[16] =  225; _outgoing1[16] =  213; _outgoing2[16] = -213; 
  _coupling[16] = 16.5/GeV; _maxweight[16] = 26.; 
  // K_2-> K* rho
  _incoming[17] =  315; _outgoing1[17] =  113; _outgoing2[17] =  313; 
  _coupling[17] = 13.42/GeV; _maxweight[17] = 30.; 
  _incoming[18] =  315; _outgoing1[18] = -213; _outgoing2[18] =  323; 
  _coupling[18] = 18.98/GeV; _maxweight[18] = 30.; 
  _incoming[19] =  325; _outgoing1[19] =  113; _outgoing2[19] =  323; 
  _coupling[19] = 13.42/GeV; _maxweight[19] = 30.; 
  _incoming[20] =  325; _outgoing1[20] =  213; _outgoing2[20] =  313; 
  _coupling[20] = 18.98/GeV; _maxweight[20] = 30.; 
  _incoming[21] = 445; _outgoing1[21] = 223; _outgoing2[21] = 223; 
  _coupling[21] = 0.00389/GeV; _maxweight[21] = 12.; 
  // initial size of the vectors
  _initsize = _incoming.size();
  // intermediates
  generateIntermediates(false);
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


void TensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << _incoming << _outgoing1 << _outgoing2 << _maxweight << ounit(_coupling,1/GeV);
}

void TensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> iunit(_coupling,1/GeV);
}

ClassDescription<TensorMesonVectorVectorDecayer> 
TensorMesonVectorVectorDecayer::initTensorMesonVectorVectorDecayer;
// Definition of the static class description member.

void TensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorVectorDecayer> documentation
    ("The TensorMesonVectorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &TensorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &TensorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &TensorMesonVectorVectorDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 1000./GeV, false, false, true);

  static ParVector<TensorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &TensorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// matrix elememt for the process
double TensorMesonVectorVectorDecayer::me2(const int,const Particle & inpart,
					   const ParticleVector & decay,
					   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  // photons ??
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix]=decay[ix]->id()==ParticleID::gamma;
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
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix],
					    outgoing,true,photon[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::calculateWaveFunctions(_vectors[ix],decay[ix],
					       outgoing,photon[ix]);
  // compute some useful dot products etc
  complex<Energy> p1eps2[3],p2eps1[3];
  Energy2 p1p2(decay[0]->momentum()*decay[1]->momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    p1eps2[ix]=_vectors[1][ix]*decay[0]->momentum();
    p2eps1[ix]=_vectors[0][ix]*decay[1]->momentum();
  }
  // compute the traces and useful dot products
  Complex trace[5];
  complex<Energy2> pboth[5];
  LorentzPolarizationVectorE pleft[2][5],pright[2][5];
  for(unsigned int ix=0;ix<5;++ix) {
    trace[ix]=_tensors[ix].xx() + _tensors[ix].yy() + _tensors[ix].zz() + _tensors[ix].tt();
    pleft[0][ix] =(-decay[0]->momentum())*_tensors[ix];
    pleft[1][ix] =(-decay[1]->momentum())*_tensors[ix];
    pright[0][ix]=_tensors[ix]*(-decay[0]->momentum());
    pright[1][ix]=_tensors[ix]*(-decay[1]->momentum());
    pboth[ix]=((pleft[0][ix]+pright[0][ix])*decay[1]->momentum());
  }
  // loop to compute the matrix element
  Complex e1e2;
  LorentzTensor<Energy2> te1e2;
  InvEnergy2 fact(_coupling[imode()]/inpart.mass());
  complex<Energy2> me;
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned iy=0;iy<3;++iy) {
      e1e2=_vectors[0][ix].dot(_vectors[1][iy]);
      te1e2=complex<Energy2>(p1p2)*
	(LorentzTensor<double>(_vectors[0][ix],_vectors[1][iy])+
	 LorentzTensor<double>(_vectors[1][iy],_vectors[0][ix]));
      for(unsigned int inhel=0;inhel<5;++inhel) {
	me = (_tensors[inhel]*te1e2
	      -p2eps1[ix]*(_vectors[1][iy].dot(pleft[0][inhel]+pright[0][inhel])) 
	      -p1eps2[iy]*(_vectors[0][ix].dot(pleft[1][inhel]+pright[1][inhel]))
	      +pboth[inhel]*e1e2
	      +(p2eps1[ix]*p1eps2[iy]-e1e2*p1p2)*trace[inhel]);
	(*ME())(inhel,ix,iy)=Complex(fact*me);
      }    
    }
  }
//   // testing the matrix element
//   double metest = ME()->contract(_rho).real();
//   Energy2 m02(sqr(inpart.mass())),m12(sqr(decay[0]->mass())),m22(sqr(decay[1]->mass()));
//   Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					     decay[1]->mass());
//   Energy2 pcm2(sqr(pcm));
//   double test = 4./15./m02/m02*sqr(_coupling[imode()])*
//     (3.*m02*(8.*pcm2*pcm2+5.*(m12*m22+pcm2*(m12+m22)))-5.*(m12-m22)*(m12-m22)*pcm2);
//   cout << "testing the matrix element VV " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << metest << " " << test <<  " " << (metest-test)/(metest+test) << endl;
  // return the answer
  return ME()->contract(_rho).real();
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
    if(id   ==_incoming[ix]) {
      if(id1==_outgoing1[ix]&&id2==_outgoing2[ix]) {
	imode=ix;
      }
      if(id2==_outgoing1[ix]&&id1==_outgoing2[ix]) {
	imode=ix;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoing1[ix]&&id2bar==_outgoing2[ix]) {
	imode=ix;
      }
      if(id2bar==_outgoing1[ix]&&id1bar==_outgoing2[ix]) {
	imode=ix;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling=_coupling[imode]*dm.parent()->mass();
  mecode=9;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}

void TensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
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
