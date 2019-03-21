// -*- C++ -*-
//
// ScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarVectorVectorDecayer class.
//

#include "ScalarVectorVectorDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix] = mode(ix)->maxWeight();
  }
}

ScalarVectorVectorDecayer::ScalarVectorVectorDecayer() 
  : _incoming(13), _outgoing1(13), _outgoing2(13), _coupling(13), 
    _maxweight(13) {
  // f_0(1370) to rho rho
  _incoming[0] = 10221; _outgoing1[0] = 113; _outgoing2[0] = 113; 
  _coupling[0] = 11.26/GeV; _maxweight[0] = 20.; 
  _incoming[1] = 10221; _outgoing1[1] = 213; _outgoing2[1] = -213; 
  _coupling[1] = 15.92/GeV; _maxweight[1] = 20.; 
  // f_0(1500) to rho rho
  _incoming[2] = 9030221; _outgoing1[2] = 113; _outgoing2[2] = 113; 
  _coupling[2] = 1.691/GeV; _maxweight[2] = 20.; 
  _incoming[3] = 9030221; _outgoing1[3] = 213; _outgoing2[3] = -213; 
  _coupling[3] = 2.391/GeV; _maxweight[3] = 20.; 
  // chi_c0 decays
  _incoming[4] = 10441; _outgoing1[4] = 443; _outgoing2[4] = 22; 
  _coupling[4] = 0.251/GeV; _maxweight[4] = 1.; 
  _incoming[5] = 10441; _outgoing1[5] = 323; _outgoing2[5] = -323; 
  _coupling[5] = 0.0088/GeV; _maxweight[5] = 1.; 
  _incoming[6] = 10441; _outgoing1[6] = 313; _outgoing2[6] = -313; 
  _coupling[6] = 0.0088/GeV; _maxweight[6] = 1.; 
  _incoming[7] = 10441; _outgoing1[7] = 333; _outgoing2[7] = 333; 
  _coupling[7] = 0.0067/GeV; _maxweight[7] = 1.; 
  _incoming[8] = 10441; _outgoing1[8] = 22; _outgoing2[8] = 22; 
  _coupling[8] = 0.0027/GeV; _maxweight[8] = 1.; 
  _incoming[12] = 10441; _outgoing1[12] = 223; _outgoing2[12] = 223; 
  _coupling[12] = 0.0093/GeV; _maxweight[12] = 1.; 
  // a'_0 -> omega rho
  _incoming[9] = 10111; _outgoing1[9] = 113; _outgoing2[9] = 223; 
  _coupling[9] = 27.09/GeV; _maxweight[9] = 20.;
  _incoming[10] = 10211; _outgoing1[10] = 213; _outgoing2[10] = 223; 
  _coupling[10] = 27.09/GeV; _maxweight[10] = 20.;
  _incoming[11] =-10211; _outgoing1[11] =-213; _outgoing2[11] = 223; 
  _coupling[11] = 27.09/GeV; _maxweight[11] = 20.; 
  // size of arrays
  _initsize = _incoming.size();
  // intermediates
  generateIntermediates(false);
}

void ScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(_coupling.size());
  if(isize!=_incoming.size()  || isize!=_outgoing1.size()||
     isize!=_outgoing2.size() || isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "ScalarVectorVectorDecayerDecayer" 
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

int ScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  cc = false;
  // check that at least some modes exist
  // must be two outgoing particles
  if(_incoming.size()==0||children.size()!=2) return -1;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  int imode(-1);
  do {
    if(_incoming[ix]==id0) {
      if((_outgoing1[ix]==id1&&_outgoing2[ix]==id2)||
	 (_outgoing1[ix]==id2&&_outgoing2[ix]==id1)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<_incoming.size());
  return imode;
}

void ScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coupling,1/GeV) << _incoming << _outgoing1 << _outgoing2 << _maxweight;
}

void ScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coupling,1/GeV) >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight;
}

ClassDescription<ScalarVectorVectorDecayer> 
ScalarVectorVectorDecayer::initScalarVectorVectorDecayer;
// Definition of the static class description member.

void ScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<ScalarVectorVectorDecayer> documentation
    ("The ScalarVectorVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static ParVector<ScalarVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &ScalarVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,int> interfaceOutcoming1
    ("FirstOutgoing",
     "The PDG code for the first outgoing particle",
     &ScalarVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,int> interfaceOutcoming2
    ("SecondOutgoing",
     "The PDG code for the second outgoing particle",
     &ScalarVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &ScalarVectorVectorDecayer::_coupling,
     1/GeV, 0, ZERO, ZERO, 10000/GeV, false, false, true);

  static ParVector<ScalarVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &ScalarVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 500000., false, false, true);

}

double ScalarVectorVectorDecayer::me2(const int,
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
  Energy2 p1p2(decay[0]->momentum()*decay[1]->momentum());
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix) {
    for(iy=0;iy<3;++iy) {
      (*ME())(0,ix,iy)=Complex(fact*(p1p2*_vectors[0][ix].dot(_vectors[1][iy])-
				     (_vectors[1][iy]*decay[0]->momentum())*
				     (_vectors[0][ix]*decay[1]->momentum())));
    }
  }
  // test of the matrix element
  //   double me = newME.contract(rhoin).real();
  //   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
  // 					   decay[1]->mass());
  //   double test = sqr(_coupling[imode()]/inpart.mass())*
  //     (2.*sqr(pcm*inpart.mass())+3.*sqr(decay[0]->mass()*decay[1]->mass()));
  //   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
  //        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
  //        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  return ME()->contract(_rho).real();
}

// output the setup info for the particle database
void ScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
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

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  itype = 12;
  return id1==_outgoing1[imode]&&id2==_outgoing2[imode];
}
