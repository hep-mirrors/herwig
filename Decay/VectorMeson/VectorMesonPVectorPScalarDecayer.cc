// -*- C++ -*-
//
// VectorMesonPVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPVectorPScalarDecayer class.
//

#include "VectorMesonPVectorPScalarDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonPVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

void VectorMesonPVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_outgoingA.size()||isize!=_outgoingP.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMesonPVectorPScalarDecayer::doinit()" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  tPDVector extpart(3);
  DecayPhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0]=getParticleData(_incoming[ix]);
    extpart[1]=getParticleData(_outgoingA[ix]);
    extpart[2]=getParticleData(_outgoingP[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

VectorMesonPVectorPScalarDecayer::VectorMesonPVectorPScalarDecayer()
  :  _coupling(21), _incoming(21), _outgoingA(21), _outgoingP(21), 
     _maxweight(21) {
  // Jpsi to K_1 K
  _incoming[0] = 443; _outgoingA[0] =  20313; _outgoingP[0] = -311; 
  _coupling[0] = 0.00127/GeV; _maxweight[0] = 12.; 
  _incoming[1] = 443; _outgoingA[1] =  20323; _outgoingP[1] = -321; 
  _coupling[1] = 0.00127/GeV; _maxweight[1] = 12.; 
  // Jpsi to b_1 pi
  _incoming[2] = 443; _outgoingA[2] =  10213; _outgoingP[2] = -211; 
  _coupling[2] = 0.00106/GeV; _maxweight[2] = 10.5; 
  _incoming[3] = 443; _outgoingA[3] =  10113; _outgoingP[3] =  111; 
  _coupling[3] = 0.00106/GeV; _maxweight[3] = 10.5; 
  // psi(2s) to K_1 K
  _incoming[4] = 100443; _outgoingA[4] =  10313; _outgoingP[4] = -311; 
  _coupling[4] = 0.00152/GeV; _maxweight[4] = 12.; 
  _incoming[5] = 100443; _outgoingA[5] =  10323; _outgoingP[5] = -321; 
  _coupling[5] = 0.00152/GeV; _maxweight[5] = 12.; 
  // psi(2s) to b_1 pi
  _incoming[6] = 100443; _outgoingA[6] =  10213; _outgoingP[6] = -211; 
  _coupling[6] = 0.000694/GeV; _maxweight[6] = 10.5; 
  _incoming[7] = 100443; _outgoingA[7] =  10113; _outgoingP[7] =  111; 
  _coupling[7] = 0.000694/GeV; _maxweight[7] = 10.5; 
  // rho'' decays
  // to h_1
  _incoming[8] =  30213; _outgoingA[8] =  10223; _outgoingP[8] = 211; 
  _coupling[8] = 1.45/GeV; _maxweight[8] = 5.5; 
  _incoming[9] =  30113; _outgoingA[9] =  10223; _outgoingP[9] = 111; 
  _coupling[9] = 1.45/GeV; _maxweight[9] = 5.5; 
  // to a_1
  _incoming[10] =  30213; _outgoingA[10] =  20213; _outgoingP[10] =  111; 
  _coupling[10] = 1.09/GeV; _maxweight[10] = 4.; 
  _incoming[11] =  30213; _outgoingA[11] =  20113; _outgoingP[11] =  211; 
  _coupling[11] = 1.09/GeV; _maxweight[11] = 4.; 
  _incoming[12] =  30113; _outgoingA[12] =  20213; _outgoingP[12] = -211; 
  _coupling[12] = 1.09/GeV; _maxweight[12] = 4.; 
  //  rho' decays
  // to h_1
  _incoming[13] =  100213; _outgoingA[13] =  10223; _outgoingP[13] = 211; 
  _coupling[13] = 1.20/GeV; _maxweight[13] = 5.; 
  _incoming[14] =  100113; _outgoingA[14] =  10223; _outgoingP[14] = 111; 
  _coupling[14] = 1.20/GeV; _maxweight[14] = 5.; 
  // to a_1
  _incoming[15] =  100213; _outgoingA[15] =  20213; _outgoingP[15] = 111; 
  _coupling[15] = 1.83/GeV; _maxweight[15] = 4.; 
  _incoming[16] =  100213; _outgoingA[16] =  20113; _outgoingP[16] = 211; 
  _coupling[16] = 1.83/GeV; _maxweight[16] = 4.; 
  _incoming[17] =  100113; _outgoingA[17] =  20213; _outgoingP[17] = -211; 
  _coupling[17] = 1.83/GeV; _maxweight[17] = 4.; 
  // omega' to b pi
  _incoming[18] = 100223; _outgoingA[18] =  10113; _outgoingP[18] =  111; 
  _coupling[18] = 1.659/GeV; _maxweight[18] = 7.; 
  _incoming[19] = 100223; _outgoingA[19] =  10213; _outgoingP[19] = -211; 
  _coupling[19] = 1.659/GeV; _maxweight[19] = 6.; 
  // psi(2s) -> h_c pi0
  _incoming[20] = 100443; _outgoingA[20] =  10443; _outgoingP[20] = 111; 
  _coupling[20] = 0.0029/GeV; _maxweight[20] = 6.; 
  // initial size of the arrays
  _initsize = _coupling.size();
  // intermediates
  generateIntermediates(false);
}

int VectorMesonPVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
      if((id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix])||
	 (id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix])) imode=ix;
    }
    if(idbar==_incoming[ix]) {
      if((id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix])||
	 (id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix])) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  return imode;
}


void VectorMesonPVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoingA << _outgoingP << _maxweight << ounit(_coupling,1/GeV);
}

void VectorMesonPVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoingA >> _outgoingP >> _maxweight >> iunit(_coupling,1/GeV);
}

ClassDescription<VectorMesonPVectorPScalarDecayer> 
VectorMesonPVectorPScalarDecayer::initVectorMesonPVectorPScalarDecayer;
// Definition of the static class description member.

void VectorMesonPVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonPVectorPScalarDecayer> documentation
    ("The VectorMesonPVectorPScalarDecayer class is designed for the "
     "decay of a vector meson to a pseudovector meson and a "
     "pseudoscalar meson.");

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonPVectorPScalarDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceOutcomingVector
    ("OutgoingPVector",
     "The PDG code for the outgoing spin-1 particle",
     &VectorMesonPVectorPScalarDecayer::_outgoingA,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,int> interfaceOutcomingScalar
    ("OutgoingPScalar",
     "The PDG code for the outgoing spin-0 particle",
     &VectorMesonPVectorPScalarDecayer::_outgoingP,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonPVectorPScalarDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 100./GeV, false, false, true);

  static ParVector<VectorMesonPVectorPScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonPVectorPScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);
}

double VectorMesonPVectorPScalarDecayer::me2(const int,
					     const Particle & inpart,
					     const ParticleVector & decay,
					     MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  // is the vector massless
  bool photon(_outgoingA[imode()]==ParticleID::gamma);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors[0],_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    VectorWaveFunction::constructSpinInfo(_vectors[1],decay[0],
					  outgoing,true,photon);
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
    return 0.;
  }
  VectorWaveFunction::calculateWaveFunctions(_vectors[1],decay[0],outgoing,photon);
  // compute the matrix element
  Energy2 p0dotpv(inpart.momentum()*decay[0]->momentum());
  complex<Energy> epsdot(ZERO);
  InvEnergy2 pre(_coupling[imode()]/inpart.mass());  
  for(unsigned int ix=0;ix<3;++ix) {
    epsdot=_vectors[1][ix]*inpart.momentum();
    for(unsigned int iy=0;iy<3;++iy) {
      (*ME())(iy,ix,0)=Complex(pre*(p0dotpv*(_vectors[1][ix].dot(_vectors[0][iy]))-
				    epsdot*(_vectors[0][iy]*decay[0]->momentum())));
    }
  }
  // test of the matrix element
//   double me = ME()->contract(_rho).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = sqr(_coupling[imode()])/3.*(2.*sqr(pcm)+3.*sqr(decay[0]->mass()));
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return ME()->contract(_rho).real();
}

bool VectorMesonPVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
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
    if(id==_incoming[ix]) {
      if(id1   ==_outgoingA[ix]&&id2   ==_outgoingP[ix]) {
	imode=ix;
	order=true;
      }
      if(id2   ==_outgoingA[ix]&&id1   ==_outgoingP[ix]) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==_incoming[ix]&&imode<0) {
      if(id1bar==_outgoingA[ix]&&id2bar==_outgoingP[ix]) {
	imode=ix;
	order=true;
      }
      if(id2bar==_outgoingA[ix]&&id1bar==_outgoingP[ix]) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<_incoming.size()&&imode<0);
  coupling = _coupling[imode]*dm.parent()->mass();  
  mecode = 4;
  return order;
}

void VectorMesonPVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " "
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":OutgoingPVector " << ix << " "
	     << _outgoingA[ix] << "\n";
      output << "newdef " << name() << ":OutgoingPScalar " << ix << " "
	     << _outgoingP[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " "
	     << _coupling[ix]*GeV << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming "  << ix << " "
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":OutgoingPVector " << ix << " "
	     << _outgoingA[ix] << "\n";
      output << "insert " << name() << ":OutgoingPScalar " << ix << " "
	     << _outgoingP[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " "
	     << _coupling[ix]*GeV << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " "
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
