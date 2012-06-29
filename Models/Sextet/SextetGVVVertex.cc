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

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
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

void SextetGVVVertex::setCoupling(Energy2 q2, tcPDPtr , tcPDPtr , 
				  tcPDPtr ) {
  if(q2 != q2Last_ || coupLast_ == 0.) {
    q2Last_ = q2;
    coupLast_ = strongCoupling(q2);
  }
  norm(coupLast_);
}
