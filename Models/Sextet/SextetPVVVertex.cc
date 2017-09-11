// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetPVVVertex class.
//

#include "SextetModel.h"
#include "SextetPVVVertex.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr SextetPVVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetPVVVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetPVVVertex,Helicity::VVVVertex>
describeSextetPVVVertex("Herwig::SextetPVVVertex", "HwSextetModel.so");

void SextetPVVVertex::Init() {

  static ClassDocumentation<SextetPVVVertex> documentation
    ("The SextetPVVVertex class implements the coupling of the gluon to two"
     " vector sextet particles");

}

void SextetPVVVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetPVVVertex::doinit()"
			       << Exception::runerror;
  if(model->VectorDoubletY16Enabled()) {
    addToList(22,ParticleID::VectorDQY16P,
	         ParticleID::VectorDQY16Pbar);
    addToList(22,ParticleID::VectorDQY16M,
	         ParticleID::VectorDQY16Mbar);

  }
  if(model->VectorDoubletY56Enabled()) {
    addToList(22,ParticleID::VectorDQY56P,
	         ParticleID::VectorDQY56Pbar);
    addToList(22,ParticleID::VectorDQY56M,
	         ParticleID::VectorDQY56Mbar);
  }
  VVVVertex::doinit();
}

void SextetPVVVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2, 
				  tcPDPtr p3) {
  if(q2 != q2Last_ || coupLast_ == 0.) {
    q2Last_ = q2;
    coupLast_ = electroMagneticCoupling(q2);
  }
  if(p1->id()==ParticleID::gamma) {
    norm(p3->iCharge()/3.*coupLast_);
  }
  else if(p2->id()==ParticleID::gamma) {
    norm(p1->iCharge()/3.*coupLast_);
  }
  else if(p3->id()==ParticleID::gamma) {
    norm(p2->iCharge()/3.*coupLast_);
  }
  else
    assert(false);
}
