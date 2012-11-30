// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFGVertex class.
//

#include "LHFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<LHFFGVertex,FFVVertex>
describeHerwigLHFFGVertex("Herwig::LHFFGVertex", "HwLHModel.so");

void LHFFGVertex::Init() {

  static ClassDocumentation<LHFFGVertex> documentation
    ("The LHFFGVertex class implements the coupling of the quarks"
     " to the gluon in the Little Higgs model");

}

// coupling for FFG vertex
void LHFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // check allowed
  int iferm=abs(a->id());
  assert((iferm>=1 && iferm<=6) || iferm==8);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  left(1.);
  right(1.);
}

void LHFFGVertex::doinit() {
  // PDG codes for the particles
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 21);
  }
  addToList(-8, 8, 21);
  FFVVertex::doinit();
}

LHFFGVertex::LHFFGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
}
