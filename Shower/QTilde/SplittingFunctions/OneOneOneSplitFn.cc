// -*- C++ -*-
//
// OneOneOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneOneOneSplitFn class.
//

#include "OneOneOneSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<OneOneOneSplitFn,Herwig::Sudakov1to2FormFactor>
describeOneOneOneSplitFn ("Herwig::OneOneOneSplitFn","HwShower.so");

void OneOneOneSplitFn::Init() {

  static ClassDocumentation<OneOneOneSplitFn> documentation
    ("The OneOneOneSplitFn class implements the g -> gg splitting function");

}
