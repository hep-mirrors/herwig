// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFGVertex class.
//

#include "LHTPFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPFFGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFGVertex::fullclone() const {
  return new_ptr(*this);
}

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<LHTPFFGVertex,FFVVertex>
describeHerwigLHTPFFGVertex("Herwig::LHTPFFGVertex", "HwLHTPModel.so");

void LHTPFFGVertex::Init() {

  static ClassDocumentation<LHTPFFGVertex> documentation
    ("The LHTPFFGVertex class implements the couples of the fermions "
     "to the gluons in the Little Higgs model with T-parity.");

}

LHTPFFGVertex::LHTPFFGVertex() 
  : coupLast_(0.), q2Last_(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

void LHTPFFGVertex::doinit() {
  // SM quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 21);
  }
  // additional top quark
  addToList(-8,  8, 21);
  // T odd quarks
  for(long ix = 4000001; ix <= 4000006; ++ix) {
    addToList(-ix,     ix, 21);
  }
  addToList(-4000008,   4000008, 21);
  FFVVertex::doinit();
}

// coupling for FFG vertex
void LHTPFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -strongCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  int iferm=abs(a->id());
  if( iferm > 8 ) iferm -= 4000000;
  if((iferm>=1 && iferm<=8)) {
    left (1.);
    right(1.);
  }
  else
    assert(false);
}
