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

MergerBase::~MergerBase() {}

void MergerBase::Init() {

  static ClassDocumentation<MergerBase> documentation
    ("MergerBase is the base class for merging helpers.");

}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractNoPIOClass<MergerBase,HandlerBase>
describeHerwigMergerBase("Herwig::MergerBase", "Herwig.so");

