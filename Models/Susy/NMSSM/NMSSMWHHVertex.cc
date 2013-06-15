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

NMSSMWHHVertex::NMSSMWHHVertex() : _sinb(0.), _cosb(0.), _sw(0.), _cw(0.),
				   _q2last(0.*MeV2), _couplast(0.) {
  orderInGem(1);
  orderInGs(0);
}

void NMSSMWHHVertex::doinit() {
  // codes for the neutral higgs
  //CP even
  int ieven[3]={25,35,45};
  //CP odd
  int iodd [2]={36,46};
  // Z CP even CP odd
  for(unsigned int ix=0;ix<3;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      addToList( 23, ieven[ix], iodd[iy] );

  // W H+ CP even
  for(unsigned int ix=0;ix<3;++ix)
    addToList( -24, 37, ieven[ix] );

   // W+ H- CP even
  for(unsigned int ix=0;ix<3;++ix)
    addToList( 24, -37, ieven[ix] );

  // W H+ CP odd
  for(unsigned int ix=0;ix<2;++ix)
    addToList( -24, 37, iodd[ix] );

  //W+ H- CP odd
  for(unsigned int ix=0;ix<2;++ix)
    addToList( 24, -37, iodd[ix] );

  // Charged higgs Z/gamma
  addToList( 22, 37, -37 );
  addToList( 23, 37, -37 );
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  // sin theta_W
  double sw2 = sin2ThetaW();
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
  _sinb = sin(beta);
  _cosb = cos(beta);
  // base class
  VSSVertex::doinit();
}

void NMSSMWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << _sinb << _cosb << _sw << _cw << _mixS << _mixP;
}

void NMSSMWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sinb >> _cosb >> _sw >> _cw >> _mixS >> _mixP;
}

ClassDescription<NMSSMWHHVertex> NMSSMWHHVertex::initNMSSMWHHVertex;
// Definition of the static class description member.

void NMSSMWHHVertex::Init() {

  static ClassDocumentation<NMSSMWHHVertex> documentation
    ("The NMSSMWHHVertex class implements the coupling of an electroweak"
     " gauge boson with two Higgs bosons in the NMSSM.");

}

//calulate the couplings
void NMSSMWHHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // weak coupling
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  // gauge bosons
  int ibos= a->id();
  int ih1 = b->id();
  int ih2 = c->id();
  Complex fact;
  if(ibos==ParticleID::Z0) {
    fact = 0.5/_cw;
    // Z H+ H-
    if(abs(ih1)==37) {
      fact *= (sqr(_cw)-sqr(_sw));
      if(ih1<0) fact *=-1.;  
    }
    // Z CP even CP odd
    else {
      if(ih1%10==6) {
	fact *= -1.; 
	swap(ih1,ih2);
      }
      int is = (ih1-25)/10;
      int ip = (ih2-36)/10;
      fact *= Complex(0.,1.)*((*_mixS)(is,1)*(*_mixP)(ip,1)-
			      (*_mixS)(is,0)*(*_mixP)(ip,0));
    }
  }
  // gamma CP even CP odd
  else if(ibos==ParticleID::gamma) {
    fact = ih1>0 ? _sw : -_sw;  
  }
  // W boson
  else {
    fact = 0.5; 
    if(abs(ih2)==37) {
      swap(ih1,ih2);
      fact*=-1; 
    }
    // H+ CP even
    if(ih2%5==0) {
      int is = (ih2-25)/10;
      fact *= (_cosb*(*_mixS)(is,1)-_sinb*(*_mixS)(is,0));
      if(ibos<0) fact*=-1; 
    }
    // H+ CP odd
    else {
      int ip = (ih2-36)/10;
      fact *=-Complex(0.,1.)*(_cosb*(*_mixP)(ip,1)+_sinb*(*_mixP)(ip,0));
    }
  }
  //output the coupling
  norm(_couplast*fact);
}
