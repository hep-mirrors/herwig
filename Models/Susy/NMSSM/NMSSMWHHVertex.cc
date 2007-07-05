// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWHHVertex class.
//

#include "NMSSMWHHVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMWHHVertex::NMSSMWHHVertex() {
  // PDG codes for the particles
  vector<int> first,second,third;
  // codes for the neutral higgs
  int ieven[3]={25,35,45};
  int iodd [2]={36,46};
  // Z S P
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      first.push_back(23);
      second.push_back(ieven[ix]);
      second.push_back(iodd [iy]);
    }
  }
  // W H+ S
  for(unsigned int ix=0;ix<3;++ix) {
    first.push_back(-24);
    second.push_back(37);
    third.push_back(ieven[ix]);
  }
  // W H+ P
  for(unsigned int ix=0;ix<2;++ix) {
    first.push_back(-24);
    second.push_back(37);
    third.push_back(iodd[ix]);
  }
  // charged higgs Z and photon
  first.push_back(22);
  second.push_back(37);
  third.push_back(-37);
  first.push_back(23);
  second.push_back(37);
  third.push_back(-37);
  // add list
  setList(first,second,third);
  _sinb=0.;
  _cosb=0.;
  _sw=0.;
  _cw=0.;
  _q2last=0.*MeV2;
  _couplast=0.;
}

void NMSSMWHHVertex::doinit() throw(InitException) {
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  _theSM = model;
  // sin theta_W
  double sw2=_theSM->sin2ThetaW();
  _sw = sqrt(sw2);
  _cw = sqrt(1.-sw2);
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWHHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMWHHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb=sin(beta);
  _cosb=cos(beta);
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  // base class
  VSSVertex::doinit();
}

void NMSSMWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << _sinb << _cosb << _sw << _cw << _mixS << _mixP << _theSM;
}

void NMSSMWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sinb >> _cosb >> _sw >> _cw >> _mixS >> _mixP >> _theSM;
}

ClassDescription<NMSSMWHHVertex> NMSSMWHHVertex::initNMSSMWHHVertex;
// Definition of the static class description member.

void NMSSMWHHVertex::Init() {

  static ClassDocumentation<NMSSMWHHVertex> documentation
    ("The NMSSMWHHVertex class implements the coupling of an electroweak"
     " gauge boson with two Higgs bosons in the NMSSM.");

}

void NMSSMWHHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // em coupling
  if(q2!=_q2last) {
    double alpha = _theSM->alphaEM(q2);
    _couplast = sqrt(4.0*Constants::pi*alpha);
    _q2last=q2;
  }
  // gauge bosons
  int ibos=a->id();
  int ih1 =b->id();
  int ih2 =c->id();
  Complex fact;
  if(ibos==ParticleID::Z0) {
    fact = 0.5/_sw/_cw;
    // Z H+ H-
    if(abs(ih1)==37) {
      fact = 0.5*(sqr(_cw)-sqr(_sw));
      if(ih1<0) fact *=-1.;
    }
    // Z S P
    else {
      if(ih1%10==6) {
	fact *=-1.;
	swap(ih1,ih2);
      }
      int is = (ih1-25)/10;
      int ip = (ih2-36)/10;
      fact *= Complex(0.,1.)*
	(*_mixS)(is,0)*(*_mixP)(ip,0)-(*_mixS)(is,1)*(*_mixP)(ip,1);
    }
  }
  else if(ibos==ParticleID::gamma) {
    fact = ih1>0 ? 1. : -1;
  }
  else {
    fact = 0.5/_sw; 
    if(abs(ih2)==37) {
      fact *=-1.;
      swap(ih1,ih2);
    }
    if(ibos>0) fact*=-1;
    // H+ S
    if(ih2%5==0) {
      int is = (ih1-25)/10;
      fact *= -              (_sinb*(*_mixS)(is,0)-_cosb*(*_mixS)(is,1));
    }
    // H+ P
    else {
      int ip = (ih1-36)/10;
      fact *= Complex(0.,1.)*(_sinb*(*_mixP)(ip,0)+_cosb*(*_mixP)(ip,1));
    }
  }
  setNorm(_couplast*fact);
}
