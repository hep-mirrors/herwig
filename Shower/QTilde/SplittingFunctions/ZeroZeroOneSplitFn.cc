// -*- C++ -*-
//
// PhitoPhiGSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ZeroZeroOneSplitFn class.
//

#include "ZeroZeroOneSplitFn.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeNoPIOClass<ZeroZeroOneSplitFn,Herwig::SudakovFormFactor>
describeZeroZeroOneSplitFn ("Herwig::ZeroZeroOneSplitFn","HwShower.so");

void ZeroZeroOneSplitFn::Init() {

  static ClassDocumentation<ZeroZeroOneSplitFn> documentation
    ("The ZeroZeroOneSplitFn class implements the splitting function for the "
     "radiation of a gluon by a scalar coloured particle");

}
