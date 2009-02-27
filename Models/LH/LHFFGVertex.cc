// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFGVertex class.
//

#include "LHFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

NoPIOClassDescription<LHFFGVertex> 
LHFFGVertex::initLHFFGVertex;
// Definition of the static class description member.

void LHFFGVertex::Init() {

  static ClassDocumentation<LHFFGVertex> documentation
    ("The LHFFGVertex class implements the coupling of the quarks"
     " to the gluon in the Little Higgs model");

}

// coupling for FFG vertex
void LHFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  setNorm(_couplast);
  // the left and right couplings
  int iferm=abs(a->id());
  if((iferm>=1 && iferm<=6) || iferm==8) {
    setLeft(1.);
    setRight(1.);
  }
  else
    throw HelicityConsistencyError() << "LHFFGVertex::setCoupling" 
				     << "Unknown particle in gluon vertex" 
				     << Exception::runerror;
}

LHFFGVertex::LHFFGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  // PDG codes for the particles
  vector<long> first,second,third;
  for(int ix=1;ix<9;++ix) {
    if(ix==7) ++ix;
    first.push_back(-ix);
    second.push_back(ix);
    third.push_back(21);
  }
  setList(first,second,third);
}
  
void LHFFGVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  FFVVertex::doinit();
}
