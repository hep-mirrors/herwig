// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHHHHVertex class.
//

#include "SMHHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMHHHHVertex::SMHHHHVertex()  : couplast_(0.), q2last_(ZERO) {}

IBPtr SMHHHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMHHHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SMHHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ratio_;
}

void SMHHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> ratio_;
}

void SMHHHHVertex::doinit() {
  SSSSVertex::doinit();
  ratio_ = -0.75*sqr(getParticleData(ParticleID::h0)->mass()/
		     getParticleData(ParticleID::Wplus)->mass());
}

ClassDescription<SMHHHHVertex> SMHHHHVertex::initSMHHHHVertex;
// Definition of the static class description member.

void SMHHHHVertex::Init() {

  static ClassDocumentation<SMHHHHVertex> documentation
    ("The SMHHHHVertex class implements the quartic Higgs"
     " coupling in the Standard Model");

}

void SMHHHHVertex::setCoupling(Energy2 q2,
			       tcPDPtr part1,tcPDPtr part2,
			       tcPDPtr part3,tcPDPtr part4) {
  assert(part1->id()==ParticleID::h0 && part2->id()==ParticleID::h0 &&
	 part3->id()==ParticleID::h0 && part4->id()==ParticleID::h0 );
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = weakCoupling(q2);
    q2last_=q2;
  }
  norm(couplast_*ratio_);
}
