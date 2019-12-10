// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMDMMediatorVertex class.
//

#include "DMDMMediatorVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DMModel.h"

using namespace Herwig;

DMDMMediatorVertex::DMDMMediatorVertex() : cDMmed_(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr DMDMMediatorVertex::clone() const {
  return new_ptr(*this);
}

IBPtr DMDMMediatorVertex::fullclone() const {
  return new_ptr(*this);
}

void DMDMMediatorVertex::persistentOutput(PersistentOStream & os) const {
  os << cDMmed_;
}

void DMDMMediatorVertex::persistentInput(PersistentIStream & is, int) {
  is >> cDMmed_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMDMMediatorVertex,FFVVertex>
describeHerwigDMDMMediatorVertex("Herwig::DMDMMediatorVertex", "HwDMModel.so");

void DMDMMediatorVertex::Init() {

  static ClassDocumentation<DMDMMediatorVertex> documentation
    ("The DMDMMediatorVertex class implements the couplnig of dark matter to the mediator.");

}

void DMDMMediatorVertex::setCoupling(Energy2 ,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  int iferm=abs(aa->id());
  assert(iferm==52);
  norm(cDMmed_);
  left(1.);
  right(1.);
}

void DMDMMediatorVertex::doinit() {
  FFVVertex::doinit();
  cDMModelPtr model = dynamic_ptr_cast<cDMModelPtr>(generator()->standardModel());
  cDMmed_ = model->cDMmed();
  addToList(52, 52, 32);
}
