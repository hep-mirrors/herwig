// -*- C++ -*-
//
// HalfOneHalfSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfOneHalfSplitFn class.
//

#include "HalfOneHalfSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<HalfOneHalfSplitFn,Herwig::SudakovFormFactor>
describeHalfOneHalfSplitFn ("Herwig::HalfOneHalfSplitFn","HwShower.so");

void HalfOneHalfSplitFn::Init() {

  static ClassDocumentation<HalfOneHalfSplitFn> documentation
    ("The HalfOneHalfSplitFn class implements the splitting "
     "function for q -> g q");

}
