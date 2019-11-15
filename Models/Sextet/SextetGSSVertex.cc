// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGSSVertex class.
//

#include "SextetGSSVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetGSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetGSSVertex,Helicity::VSSVertex,false,true>
describeSextetGSSVertex("Herwig::SextetGSSVertex", "HwSextetModel.so");

void SextetGSSVertex::Init() {

  static ClassDocumentation<SextetGSSVertex> documentation
    ("The SextetGSSVertex class implements the coupling of the gluon"
     " to scalar diquarks.");

}

void SextetGSSVertex::doinit() {
  orderInGs (1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(21,ParticleID::ScalarDQSingletY43,
   	         ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(21,ParticleID::ScalarDQSingletY13,
	         ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(21,ParticleID::ScalarDQSingletY23,
	         ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(21,ParticleID::ScalarDQTripletP,
	         ParticleID::ScalarDQTripletPbar);
    addToList(21,ParticleID::ScalarDQTriplet0,
	         ParticleID::ScalarDQTriplet0bar);
    addToList(21,ParticleID::ScalarDQTripletM,
	         ParticleID::ScalarDQTripletMbar);
  }
  Helicity::VSSVertex::doinit();
}

void SextetGSSVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
				  tcPDPtr part1,
#else
				  tcPDPtr ,
#endif
				  tcPDPtr part2, tcPDPtr ) {
  assert(part1->id()==ParticleID::g);
#ifndef NDEBUG
  long idq = abs(part2->id());
#endif
  assert(idq == ParticleID::ScalarDQSingletY43 ||
	 idq == ParticleID::ScalarDQSingletY13 ||
	 idq == ParticleID::ScalarDQSingletY23 ||
	 idq == ParticleID::ScalarDQTripletP   ||
	 idq == ParticleID::ScalarDQTriplet0   ||
	 idq == ParticleID::ScalarDQTripletM);
  if(q2 != q2Last_ || coupLast_ == 0.) {
    coupLast_ = strongCoupling(q2);
    q2Last_   = q2;
  }
  if(part2->id()>0) 
    norm(-coupLast_);
  else
    norm( coupLast_);
}
