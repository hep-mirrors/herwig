// -*- C++ -*-
//
// PScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorVectorDecayer class.
//
#include "PScalarVectorVectorDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      _maxweight[ix] = mode(ix)->maxWeight();
  }
}

PScalarVectorVectorDecayer::PScalarVectorVectorDecayer() 
  : _incoming(10), _outgoing1(10), _outgoing2(10), _coupling(10), _maxweight(10) {
  // decay eta -> omega gamma
  _incoming[0] = 331; _outgoing1[0] = 223; _outgoing2[0] = 22; 
  _coupling[0] = 0.1412/GeV; _maxweight[0] = 1.2; 
  // decay pi -> gamma gamma
  _incoming[1] = 111; _outgoing1[1] = 22; _outgoing2[1] = 22; 
  _coupling[1] = 0.0178/GeV; _maxweight[1] = 1.1; 
  // decay eta -> gamma gamma
  _incoming[2] = 221; _outgoing1[2] = 22; _outgoing2[2] = 22; 
  _coupling[2] = 0.0176/GeV; _maxweight[2] = 1.1; 
  // decay eta' -> gamma gamma
  _incoming[3] = 331; _outgoing1[3] = 22; _outgoing2[3] = 22; 
  _coupling[3] = 0.0221/GeV; _maxweight[3] = 1.1; 
  // decay eta_c -> rho rho
  _incoming[4] = 441; _outgoing1[4] = 213; _outgoing2[4] = -213; 
  _coupling[4] = 0.0525/GeV; _maxweight[4] = 2.7; 
  _incoming[5] = 441; _outgoing1[5] = 113; _outgoing2[5] =  113; 
  _coupling[5] = 0.0371/GeV; _maxweight[5] = 2.7; 
  // decay eta-c -> phi phi
  _incoming[6] = 441; _outgoing1[6] = 333; _outgoing2[6] = 333; 
  _coupling[6] = 0.0267/GeV; _maxweight[6] = 9.; 
  // decay eta-c -> gamma gamma
  _incoming[7] = 441; _outgoing1[7] = 22; _outgoing2[7] = 22; 
  _coupling[7] = 0.00521/GeV; _maxweight[7] = 1.2; 
  // decay eta_c -> K* K*
  _incoming[8] = 441; _outgoing1[8] = 323; _outgoing2[8] = -323; 
  _coupling[8] = 0.0308/GeV; _maxweight[8] = 5.3; 
  _incoming[9] = 441; _outgoing1[9] = 313; _outgoing2[9] = -313; 
  _coupling[9] = 0.0308/GeV; _maxweight[9] = 5.3; 
  // initial size of the vectors
  _initsize = _incoming.size();
  // intermediates
  generateIntermediates(false);
}

void PScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoing1.size()||
     isize!=_outgoing2.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  DecayPhaseSpaceModePtr mode;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    extpart[0] = getParticleData(_incoming[ix]);
    extpart[1] = getParticleData(_outgoing1[ix]);
    extpart[2] = getParticleData(_outgoing2[ix]);
    if(extpart[0]&&extpart[1]&&extpart[2]) 
      mode=new_ptr(DecayPhaseSpaceMode(extpart,this));
    else
      mode=DecayPhaseSpaceModePtr();
    addMode(mode,_maxweight[ix],wgt);
  }
}

int PScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  cc = false;
  if(children.size()!=2) return -1;
  int id(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  int imode(-1);
  do {
    if(_incoming[ix]==id) {
      if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	 (id2==_outgoing1[ix]&&id1==_outgoing2[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV) >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<PScalarVectorVectorDecayer> 
PScalarVectorVectorDecayer::initPScalarVectorVectorDecayer;
// Definition of the static class description member.

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The PScalarVectorVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<PScalarVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &PScalarVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &PScalarVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &PScalarVectorVectorDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 10000/GeV, false, false, true);

  static ParVector<PScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &PScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 200., false, false, true);
}

double PScalarVectorVectorDecayer::me2(const int,
				       const Particle & inpart,
				       const ParticleVector & decay,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&inpart),incoming);
  }
  if(meopt==Terminate) {
    // set up the spin information for the decay products
    ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),
					  incoming,true);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix],
					    outgoing,true,photon[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vectors[ix],decay[ix],outgoing,photon[ix]);
  // now compute the matrix element
  InvEnergy2 fact(_coupling[imode()]/inpart.mass());
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix) {
    for(iy=0;iy<3;++iy) {
      (*ME())(0,ix,iy)=fact*epsilon(_vectors[0][ix],decay[1]->momentum(),
				 _vectors[1][iy])
	*decay[0]->momentum();
    }
  }
  // test of the matrix element
//   double test = 2.*sqr(fact*inpart.mass())*
//     sqr(Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),decay[1]->mass()));
//   double me = newME.contract(rhoin).real();
//   cerr << "testing the matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " " 
//        << me << " " << (me-test)/(me+test) << "\n";
  return ME()->contract(_rho).real();
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do {
    if(_incoming[ix]==id) {
      if((id1==_outgoing1[ix]&&id2==_outgoing2[ix])||
	 (id2==_outgoing1[ix]&&id1==_outgoing2[ix])) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  coupling=_coupling[imode]*dm.parent()->mass();
  itype = 3;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}

// output the setup info for the particle database
void PScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming   " << ix << " "
	     << _incoming[ix]   << "\n";
      output << "newdef " << name() << ":FirstOutgoing  " << ix << " "
	     << _outgoing1[ix]  << "\n";
      output << "newdef " << name() << ":SecondOutgoing " << ix << " "
	     << _outgoing2[ix]  << "\n";
      output << "newdef " << name() << ":Coupling   " << ix << " "
	     << _coupling[ix]*GeV   << "\n";
      output << "newdef " << name() << ":MaxWeight  " << ix << " "
	     << _maxweight[ix]  << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming   " << ix << " "
	     << _incoming[ix]   << "\n";
      output << "insert " << name() << ":FirstOutgoing  " << ix << " "
	     << _outgoing1[ix]  << "\n";
      output << "insert " << name() << ":SecondOutgoing " << ix << " "
	     << _outgoing2[ix]  << "\n";
      output << "insert " << name() << ":Coupling   " << ix << " "
	     << _coupling[ix]*GeV   << "\n";
      output << "insert " << name() << ":MaxWeight  " << ix << " "
	     << _maxweight[ix]  << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
