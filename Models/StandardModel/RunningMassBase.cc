// -*- C++ -*-
//
// RunningMassBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningMassBase class.
//

#include "RunningMassBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RunningMassBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(_theMass, GeV);
}

void RunningMassBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_theMass, GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<RunningMassBase,Interfaced>
describeHerwigRunningMassBase("Herwig::RunningMassBase", "Herwig.so");

void RunningMassBase::Init() {
 
  static ClassDocumentation<RunningMassBase> documentation
    ("The RunningMassBase class is the base class for running mass"
     "calculations");
  
}

void RunningMassBase::doinit() {
  _theMass = mass();
  Interfaced::doinit();
}
