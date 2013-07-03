// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGVVVertex class.
//

#include "SextetModel.h"
#include "SextetGVVVertex.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr SextetGVVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGVVVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetGVVVertex,Helicity::VVVVertex>
describeSextetGVVVertex("Herwig::SextetGVVVertex", "HwSextetModel.so");

void SextetGVVVertex::Init() {

  static ClassDocumentation<SextetGVVVertex> documentation
    ("The SextetGVVVertex class implements the coupling of the gluon to two"
     " vector sextet particles");

}

void SextetGVVVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGVVVertex::doinit()"
			       << Exception::runerror;
  if(model->VectorDoubletY16Enabled()) {
    addToList(21,ParticleID::VectorDQY16P,
	         ParticleID::VectorDQY16Pbar);
    addToList(21,ParticleID::VectorDQY16M,
	         ParticleID::VectorDQY16Mbar);

  }
  if(model->VectorDoubletY56Enabled()) {
    addToList(21,ParticleID::VectorDQY56P,
	         ParticleID::VectorDQY56Pbar);
    addToList(21,ParticleID::VectorDQY56M,
	         ParticleID::VectorDQY56Mbar);
  }
  VVVVertex::doinit();
}

void SextetGVVVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2, 
				  tcPDPtr p3) {
  if(q2 != q2Last_ || coupLast_ == 0.) {
    q2Last_ = q2;
    coupLast_ = strongCoupling(q2);
  }
  if((p1->id()==ParticleID::g&&p2->id()>0&&p3->id()<0)||
     (p2->id()==ParticleID::g&&p3->id()>0&&p1->id()<0)||
     (p3->id()==ParticleID::g&&p1->id()>0&&p2->id()<0))
    norm(-coupLast_);
  else
    norm( coupLast_);
}
