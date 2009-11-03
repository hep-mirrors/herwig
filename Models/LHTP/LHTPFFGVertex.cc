// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFGVertex class.
//

#include "LHTPFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPFFGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFGVertex::fullclone() const {
  return new_ptr(*this);
}

NoPIOClassDescription<LHTPFFGVertex> 
LHTPFFGVertex::initLHTPFFGVertex;
// Definition of the static class description member.

void LHTPFFGVertex::Init() {

  static ClassDocumentation<LHTPFFGVertex> documentation
    ("There is no documentation for the LHTPFFGVertex class");

}

LHTPFFGVertex::LHTPFFGVertex() 
  : _couplast(0.), _q2last(0.*GeV2) {
  // SM quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 21);
  }
  // additional top quark
  addToList(-8,  8, 21);
  // T odd quarks
  for(long ix = 4000001; ix < 4000006; ++ix) {
    addToList(-ix,     ix, 21);
  }
  addToList(-4000008,   4000008, 21);
}
  
void LHTPFFGVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  FFVVertex::doinit();
}

// coupling for FFG vertex
void LHTPFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if( iferm > 8 ) iferm -= 4000000;
  if((iferm>=1 && iferm<=8)) {
    left(1.);
    right(1.);
  }
  else
    throw HelicityConsistencyError() << "LHTPFFGVertex::setCoupling" 
				     << "Unknown particle in gluon vertex" 
				     << Exception::runerror;
}
