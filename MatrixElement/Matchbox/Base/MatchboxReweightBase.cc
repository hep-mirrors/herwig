// -*- C++ -*-
//
// MatchboxReweightBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxReweightBase class.
//

#include "MatchboxReweightBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MatchboxReweightBase::cloneDependencies(const std::string&) {}

void MatchboxReweightBase::persistentOutput(PersistentOStream &) const {}

void MatchboxReweightBase::persistentInput(PersistentIStream &, int) {}

void MatchboxReweightBase::Init() {

  static ClassDocumentation<MatchboxReweightBase> documentation
    ("MatchboxReweightBase");

}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<MatchboxReweightBase,HandlerBase>
describeMatchboxReweightBase("Herwig::MatchboxReweightBase", "Herwig.so");
