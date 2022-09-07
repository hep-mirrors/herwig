// -*- C++ -*-
//
// MergerBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MergerBase class.
//

#include "MergerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"



using namespace Herwig;

MergerBase::MergerBase() 
  : HandlerBase() {}

void MergerBase::Init() {

  static ClassDocumentation<MergerBase> documentation
    ("MergerBase is the base class for merging helpers.");

}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<MergerBase,HandlerBase>
describeHerwigMergerBase("Herwig::MergerBase", "Herwig.so");

