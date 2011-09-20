// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHHHVertex class.
//

#include "SMHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMHHHVertex::SMHHHVertex() : ratio_(ZERO), couplast_(0.), q2last_(ZERO) {
  orderInGem(1);
  orderInGs (0);
}

IBPtr SMHHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMHHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SMHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(ratio_,GeV);
}

void SMHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ratio_,GeV);
}

ClassDescription<SMHHHVertex> SMHHHVertex::initSMHHHVertex;
// Definition of the static class description member.

void SMHHHVertex::Init() {

  static ClassDocumentation<SMHHHVertex> documentation
    ("The SMHHHVertex class implements the triple Higgs"
     " coupling in the Standard Model.");

}

void SMHHHVertex::doinit() {
  addToList(25,25,25);
  SSSVertex::doinit();
  ratio_ = -1.5*sqr(getParticleData(ParticleID::h0)->mass())/
    getParticleData(ParticleID::Wplus)->mass();
}

#ifndef NDEBUG
void SMHHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3) {
#else
void SMHHHVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,tcPDPtr) {
#endif
  assert(part1->id()==ParticleID::h0 &&
	 part2->id()==ParticleID::h0 &&
	 part3->id()==ParticleID::h0 );
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = weakCoupling(q2)*ratio_*UnitRemoval::InvE;
    q2last_=q2;
  }
  norm(couplast_);
}
