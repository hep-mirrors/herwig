  // -*- C++ -*-
  //
  // CMWOneOneOneSplitFn.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the CMWOneOneOneSplitFn class.
  //

#include "CMWOneOneOneSplitFn.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"


using namespace Herwig;

DescribeClass<CMWOneOneOneSplitFn,Herwig::OneOneOneSplitFn>
describeCMWOneOneOneSplitFn ("Herwig::CMWOneOneOneSplitFn","HwShower.so");

void CMWOneOneOneSplitFn::Init() {
  
  static ClassDocumentation<CMWOneOneOneSplitFn> documentation
  ("The CMWOneOneOneSplitFn class implements the g -> gg splitting function");
  
  static Switch<CMWOneOneOneSplitFn, bool> interfaceIsIS
  ("isInititalState",
   "Switch on if this kernel is used for initial state emission.",
   &CMWOneOneOneSplitFn::isIS_, 0, false, false);
  static SwitchOption interfaceIsISYes
  (interfaceIsIS,"No","The kernel is used for final state emissions.", 0);
  static SwitchOption interfaceIsISNo
  (interfaceIsIS,"Yes","The kernel is used for final state emissions.", 1);
  
}

void CMWOneOneOneSplitFn::persistentOutput(PersistentOStream & os) const {
  os << isIS_ ;
}

void CMWOneOneOneSplitFn::persistentInput(PersistentIStream & is, int) {
  is >> isIS_ ;
}
