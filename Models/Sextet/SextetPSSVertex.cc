// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetPSSVertex class.
//

#include "SextetPSSVertex.h"
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

IBPtr SextetPSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetPSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetPSSVertex,Helicity::VSSVertex,false,true>
describeSextetPSSVertex("Herwig::SextetPSSVertex", "HwSextetModel.so");

void SextetPSSVertex::Init() {

  static ClassDocumentation<SextetPSSVertex> documentation
    ("The SextetPSSVertex class implements the coupling of the gluon"
     " to scalar diquarks.");

}

void SextetPSSVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetPSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(22,ParticleID::ScalarDQSingletY43,
   	         ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(22,ParticleID::ScalarDQSingletY13,
	         ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(22,ParticleID::ScalarDQSingletY23,
	         ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(22,ParticleID::ScalarDQTripletP,
	         ParticleID::ScalarDQTripletPbar);
    addToList(22,ParticleID::ScalarDQTriplet0,
	         ParticleID::ScalarDQTriplet0bar);
    addToList(22,ParticleID::ScalarDQTripletM,
	         ParticleID::ScalarDQTripletMbar);
  }
  Helicity::VSSVertex::doinit();
}

void SextetPSSVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
				  tcPDPtr part1,
#else
				  tcPDPtr ,
#endif
				  tcPDPtr part2, tcPDPtr ) {
  assert(part1->id()==ParticleID::gamma);
  tcPDPtr sca = part2->id()>0 ? part2 : tcPDPtr(part2->CC());
  assert(sca->id() == ParticleID::ScalarDQSingletY43 ||
	 sca->id() == ParticleID::ScalarDQSingletY13 ||
	 sca->id() == ParticleID::ScalarDQSingletY23 ||
	 sca->id() == ParticleID::ScalarDQTripletP   ||
	 sca->id() == ParticleID::ScalarDQTriplet0   ||
	 sca->id() == ParticleID::ScalarDQTripletM);
  if(q2 != q2Last_ || coupLast_ == 0.) {
    coupLast_ = electroMagneticCoupling(q2);
    q2Last_   = q2;
  }
  if(part2->id()>0) 
    norm(-sca->iCharge()/3.*coupLast_);
  else
    norm( sca->iCharge()/3.*coupLast_);
}
