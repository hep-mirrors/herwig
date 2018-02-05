// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWHHVertex class.
//

#include "SMWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMWWHHVertex::SMWWHHVertex() : ratio_(0.), couplast_(0.), q2last_(ZERO) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SMWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SMWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ratio_;
}

void SMWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> ratio_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWWHHVertex,Helicity::VVSSVertex>
describeHerwigSMWWHHVertex("Herwig::SMWWHHVertex", "Herwig.so");

void SMWWHHVertex::Init() {

  static ClassDocumentation<SMWWHHVertex> documentation
    ("The SMWWHHVertex class implements the coupling of two electroweeak"
     " gauge bosons and the Higgs boson in the Standard Model.");

}

void SMWWHHVertex::doinit() {
  addToList( 23, 23, 25, 25);
  addToList( 24,-24, 25, 25);
  VVSSVertex::doinit();
  ratio_ = 1./(1.-sin2ThetaW());
}

void SMWWHHVertex::setCoupling(Energy2 q2,
			       tcPDPtr part1,tcPDPtr,
#ifndef NDEBUG
			       tcPDPtr part3,tcPDPtr part4) {
#else
			       tcPDPtr,tcPDPtr) {
#endif
  assert(part3->id()==ParticleID::h0 && part4->id()==ParticleID::h0 );
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(weakCoupling(q2));
    q2last_=q2;
  }
  if(part1->id()==ParticleID::Z0) {
    norm(0.5*couplast_*ratio_);
  }
  else {
    norm(0.5*couplast_);
  }
}
