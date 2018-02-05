// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGGSSVertex class.
//

#include "SextetGGSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "SextetModel.h"
#include "SextetParticles.h"

using namespace Herwig;

SextetGGSSVertex::SextetGGSSVertex() : q2last_(), couplast_() {
  colourStructure(ColourStructure::SU3TT6);
}

IBPtr SextetGGSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGGSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetGGSSVertex,VVSSVertex>
describeSextetGGSSVertex("Herwig::SextetGGSSVertex", "HwSextetModel.so");

void SextetGGSSVertex::Init() {

  static ClassDocumentation<SextetGGSSVertex> documentation
    ("The SextetGGSSVertex class implements the coupling of two gluons to two"
     " scalar sextets");

}

void SextetGGSSVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				   tcPDPtr) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));
    q2last_ = q2;
  }
  norm(couplast_);
}

void SextetGGSSVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY43,
	            ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY13,
	            ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY23,
	            ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(21,21,ParticleID::ScalarDQTripletP,
	            ParticleID::ScalarDQTripletPbar);
    addToList(21,21,ParticleID::ScalarDQTriplet0,
	            ParticleID::ScalarDQTriplet0bar);
    addToList(21,21,ParticleID::ScalarDQTripletM,
	            ParticleID::ScalarDQTripletMbar);
  }
  VVSSVertex::doinit();
}
