// -*- C++ -*-
//
// SSGGSQSQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSGGSQSQVertex class.
//

#include "SSGGSQSQVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG::Helicity;
using namespace Herwig;

SSGGSQSQVertex::SSGGSQSQVertex() : q2last_(),couplast_(0.) {
  colourStructure(ColourStructure::SU3TTFUNDS);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SSGGSQSQVertex,Helicity::VVSSVertex>
describeSSGGSQSQVertex("Herwig::SSGGSQSQVertex", "HwSusy.so");

void SSGGSQSQVertex::Init() {

  static ClassDocumentation<SSGGSQSQVertex> documentation
    ("This implements the gluon-gluon-squark-squark vertex.");

}

void SSGGSQSQVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				 tcPDPtr) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));
    q2last_ = q2;
  }
  norm(couplast_);
}

void SSGGSQSQVertex::doinit() {
  //L-L squarks
  for(long ix=1000001;ix<1000007;++ix) {
    addToList(21,21,ix,-ix);
  }
  //R-R squarks
  for(long ix=2000001;ix<2000007;++ix) {
    addToList(21,21,ix,-ix);
  }
  orderInGs(2);
  orderInGem(0);
  VVSSVertex::doinit();
}
