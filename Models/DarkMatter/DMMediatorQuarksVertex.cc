// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMMediatorQuarksVertex class.
//

#include "DMMediatorQuarksVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DMModel.h"

using namespace Herwig;

DMMediatorQuarksVertex::DMMediatorQuarksVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr DMMediatorQuarksVertex::clone() const {
  return new_ptr(*this);
}

IBPtr DMMediatorQuarksVertex::fullclone() const {
  return new_ptr(*this);
}

void DMMediatorQuarksVertex::persistentOutput(PersistentOStream & os) const {
  os << cSMmed_;
}

void DMMediatorQuarksVertex::persistentInput(PersistentIStream & is, int) {
  is >> cSMmed_;
}

void DMMediatorQuarksVertex::setCoupling(Energy2 ,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  int iferm=abs(aa->id());
  assert(iferm>0 && iferm<4);
  norm(cSMmed_[iferm-1]);
  left(1.);
  right(1.);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMMediatorQuarksVertex,FFVVertex>
  describeHerwigDMMediatorQuarksVertex("Herwig::DMMediatorQuarksVertex", "HwDMModel.so");

void DMMediatorQuarksVertex::Init() {

  static ClassDocumentation<DMMediatorQuarksVertex> documentation
    ("The DMMediatorQuarksVertex class implements the coupling of the quarks to the mediator.");

}

void DMMediatorQuarksVertex::doinit() {
  FFVVertex::doinit();
  cDMModelPtr model = dynamic_ptr_cast<cDMModelPtr>(generator()->standardModel());
  cSMmed_ = model->cSMmed();
  for(int ix=1;ix<4;++ix) {
    addToList(-ix, ix, 32);
  }
}
