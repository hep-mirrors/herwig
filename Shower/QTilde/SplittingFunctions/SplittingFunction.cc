// -*- C++ -*-
//
// SplittingFunction.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SplittingFunction class.
//

#include "SplittingFunction.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeAbstractNoPIOClass<SplittingFunction,Interfaced>
describeSplittingFunction ("Herwig::SplittingFunction","");

void SplittingFunction::Init() {

  static ClassDocumentation<SplittingFunction> documentation
    ("The SplittingFunction class is the based class for 1->2 splitting functions"
     " in Herwig");
}
