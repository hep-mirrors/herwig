// -*- C++ -*-
//
// HalfHalfOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HalfHalfOneSplitFn class.
//

#include "HalfHalfOneSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<HalfHalfOneSplitFn,Herwig::Sudakov1to2FormFactor>
describeHalfHalfOneSplitFn ("Herwig::HalfHalfOneSplitFn","HwShower.so");

void HalfHalfOneSplitFn::Init() {

  static ClassDocumentation<HalfHalfOneSplitFn> documentation
    ("The HalfHalfOneSplitFn class implements the q -> qg splitting function");

}
