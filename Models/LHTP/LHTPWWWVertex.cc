// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWWWVertex class.
//

#include "LHTPWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "LHTPModel.h"

using namespace Herwig;

LHTPWWWVertex::LHTPWWWVertex() : couplast_(0.), q2last_(0.*MeV2),
				 couplings_(3 ,0.) {
  orderInGem(1);
  orderInGs(0);
  //SM interactions
  addToList(24,  -24,    22);
  addToList(24,  -24,  23);

  //LHTP
  //W_H W_H A_L
  addToList(34,  -34,    22);
  //W_H W_H Z_L
  addToList(34,  -34,  23);
  //W_H W_L A_H
  addToList(34,  -24,  32);

  addToList(-34,  24,  32); 
  //W_H W_L Z_H
  addToList(34,  -24,  33);
  addToList(-34,  24,  33);
}

void LHTPWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << couplings_;
}

void LHTPWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> couplings_;
}

ClassDescription<LHTPWWWVertex> LHTPWWWVertex::initLHTPWWWVertex;
// Definition of the static class description member.

void LHTPWWWVertex::Init() {

  static ClassDocumentation<LHTPWWWVertex> documentation
    ("There is no documentation for the LHTPWWWVertex class");

}

void LHTPWWWVertex::doinit() {
  VVVVertex::doinit();
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "LHTPWWWVertex::doinit() - Model pointer must be of LHTPModel"
      << "type, cannot continue without this."
      << Exception::abortnow;
  
  double sw(sqrt(model->sin2ThetaW()));
  double cw(sqrt(1. - model->sin2ThetaW()));

  //W W Z
  couplings_[0] = cw/sw;
  //W_L W_H A_H
  couplings_[1] = -model->sinThetaH()/sw;
  //W_L W_H Z_H
  couplings_[2] = 1./sw;

}

void LHTPWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  if(q2 != q2last_) {
    couplast_ = electroMagneticCoupling(q2);
    q2last_ = q2;
  }
  long boson(0);
  double perm(1.);
  //rearrange to keep order but with non W at [2]
  if( abs(a->id()) != 24 && abs(a->id()) != 34 ) {
    boson = a->id();
    if( c->id() < 0 ) perm = -1.;
  }
  else if( abs(b->id()) != 24 && abs(b->id()) != 34 ) {
    boson = b->id();
    if( c->id() < 0 ) perm = -1.;
  }
  else if( abs(c->id()) != 24 && abs(c->id()) != 34 ) {
    boson = c->id();
    if( b->id() < 0 ) perm = -1.;
  }
  else 
    throw Helicity::HelicityConsistencyError() 
      << "LHTPWWWVertex::setCoupling - Incorrect particles in LHTPWWWVertex. " 
      << a->id() << " " << b->id() << " " << c->id() << '\n'
      << Exception::runerror;

  if( boson == 22 )
    norm(perm*couplast_);
  else if( boson == 23 )
    norm(perm*couplings_[0]*couplast_);
  else if( boson == 32 )
    norm(perm*couplings_[1]*couplast_);
  else if( boson == 33 )
    norm(perm*couplings_[2]*couplast_);
  else 
    throw Helicity::HelicityConsistencyError() 
      << "LHTPWWWVertex::setCoupling - Incorrect boson in LHTPWWWVertex. "
      << a->id() << " " << b->id() << " " << c->id() << '\n'
      << Exception::runerror;
}
