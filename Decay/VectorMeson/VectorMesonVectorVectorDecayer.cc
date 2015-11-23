// -*- C++ -*-
//
// VectorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorVectorDecayer class.
//

#include "VectorMesonVectorVectorDecayer.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_incoming.size();++ix)
      if(mode(ix)) _maxweight[ix]=mode(ix)->maxWeight();
  }
}

void VectorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  unsigned int isize(_incoming.size());
  if(isize!=_outgoing1.size()||isize!=_outgoing2.size()||
     isize!=_maxweight.size()||isize!=_coupling.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorVectorDecayer" << Exception::runerror;
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

VectorMesonVectorVectorDecayer::VectorMesonVectorVectorDecayer() 
  : _coupling(4), _incoming(4), _outgoing1(4), _outgoing2(4), 
    _maxweight(4) {
  // decay of rho'' to rho rho
  _incoming[0] = 30213; _outgoing1[0] =  213; _outgoing2[0] = 113; 
  _coupling[0] = 3.21; _maxweight[0] = 35.; 
  _incoming[1] = 30113; _outgoing1[1] = -213; _outgoing2[1] = 213; 
  _coupling[1] = 3.21; _maxweight[1] = 22.; 
  // decay of rho' to rho rho
  _incoming[2] =  100213; _outgoing1[2] =  213; _outgoing2[2] = 113; 
  _coupling[2] = 9.59; _maxweight[2] = 55.; 
  _incoming[3] =  100113; _outgoing1[3] = -213; _outgoing2[3] = 213; 
  _coupling[3] = 9.59; _maxweight[3] = 50.; 
  // initial size of the arrays
  _initsize=_incoming.size();
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

void VectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << _incoming << _outgoing1 << _outgoing2 << _maxweight << _coupling;
}

void VectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _incoming >> _outgoing1 >> _outgoing2 >> _maxweight >> _coupling;
}

ClassDescription<VectorMesonVectorVectorDecayer> 
VectorMesonVectorVectorDecayer::initVectorMesonVectorVectorDecayer;
// Definition of the static class description member.

void VectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorVectorDecayer> documentation
    ("The VectorMesonVectorVectorDecayer class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &VectorMesonVectorVectorDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing1
    ("Outgoing1",
     "The PDG code for the first outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing1,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,int> interfaceOutgoing2
    ("Outgoing2",
     "The PDG code for the second outgoing particle",
     &VectorMesonVectorVectorDecayer::_outgoing2,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &VectorMesonVectorVectorDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<VectorMesonVectorVectorDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &VectorMesonVectorVectorDecayer::_maxweight,
     0, 0, 0, 0., 1000., false, false, true);

}

double VectorMesonVectorVectorDecayer::me2(const int,
					   const Particle & inpart,
					   const ParticleVector & decay,
					   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix) 
    photon[ix] = decay[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors[0],_rho,
					       const_ptr_cast<tPPtr>(&inpart),
					       incoming,false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::constructSpinInfo(_vectors[0],const_ptr_cast<tPPtr>(&inpart),
					  incoming,true,false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::constructSpinInfo(_vectors[ix+1],decay[ix],
					    outgoing,true,photon[ix]);
    return 0.;
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(_vectors[ix+1],decay[ix],outgoing,photon[ix]);
  // work out the dot products we need for the matrix element
  Energy2 p1p2((decay[0]->momentum())*(decay[1]->momentum()));
  complex<Energy> p1eps2[3],p2eps1[3];
  for(unsigned int ix=0;ix<3;++ix) {
    p1eps2[ix]=_vectors[2][ix]*(decay[0]->momentum());
    p2eps1[ix]=_vectors[1][ix]*(decay[1]->momentum());
  }
  // compute the matrix element
  Lorentz5Momentum pdiff(decay[0]->momentum()-decay[1]->momentum());
  Energy2 m12(decay[0]->mass()*decay[0]->mass()),m22(decay[1]->mass()*decay[1]->mass());
  InvEnergy3 fact(2.*_coupling[imode()]/(inpart.mass()*inpart.mass()*inpart.mass()));
  LorentzPolarizationVector vtemp;
  for(unsigned int ipol1=0;ipol1<3;++ipol1) {
    for(unsigned int ipol2=0;ipol2<3;++ipol2) {
      Complex eps1eps2=_vectors[1][ipol1].dot(_vectors[2][ipol2]);
      vtemp=fact*(p1eps2[ipol2]*p2eps1[ipol1]*pdiff
		  +p1eps2[ipol2]*m22*_vectors[1][ipol1]
		  -p2eps1[ipol1]*m12*_vectors[2][ipol2]
		  +eps1eps2*(-p1p2*pdiff+m12*decay[1]->momentum()
			     -m22*decay[0]->momentum()));
      for(unsigned int inpol=0;inpol<3;++inpol) 
	(*ME())(inpol,ipol1,ipol2)=_vectors[0][inpol].dot(vtemp);
    }
  }
  // test of the matrix element
//   double me = newME.contract(rhoin).real();
//   Energy pcm=Kinematics::pstarTwoBodyDecay(inpart.mass(),decay[0]->mass(),
// 					   decay[1]->mass());
//   double test = 8./3.*sqr(_coupling[imode()]*pcm/inpart.mass())*
//     (1.+sqr(decay[0]->mass()/inpart.mass())+sqr(decay[1]->mass()/inpart.mass()));
//   cerr << "testing matrix element for " << inpart.PDGName() << " -> " 
//        << decay[0]->PDGName() << " " << decay[1]->PDGName() << " "
//        << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return ME()->contract(_rho).real();
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
  coupling = _coupling[imode]; 
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
  for(unsigned int ix=0;ix<_incoming.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "newdef " << name() << ":Outgoing1 " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "newdef " << name() << ":Outgoing2 " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "newdef " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":Incoming " << ix << " " 
	     << _incoming[ix] << "\n";
      output << "insert " << name() << ":Outgoing1 " << ix << " " 
	     << _outgoing1[ix] << "\n";
      output << "insert " << name() << ":Outgoing2 " << ix << " " 
	     << _outgoing2[ix] << "\n";
      output << "insert " << name() << ":Coupling " << ix << " " 
	     << _coupling[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " << ix << " " 
	     << _maxweight[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
