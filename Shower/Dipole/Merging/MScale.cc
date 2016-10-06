  // -*- C++ -*-
  //
  // MScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2012 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the MScale class.
  //

#include "MScale.h"
#include "Node.h"


using namespace Herwig;

MScale::MScale() {}

MScale::~MScale() {}



Energy2 MScale::renormalizationScale() const {

  if (theMergingHelper->renormscale() != ZERO ){
      return sqr(theMergingHelper->renormscale());
  }
  if (!theScaleChoice){
    cerr<<"MScale: You need to set a scale for MScale.";
    assert(false);
  }
  theScaleChoice->setXComb(theLastXComb);
  return theScaleChoice->renormalizationScale();
}

Energy2 MScale::factorizationScale() const {
  if ( theMergingHelper->renormscale() != ZERO )
      return sqr( theMergingHelper->renormscale());
  
  return sqr(10.*GeV); // TODO fix in PDF workflow
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

#include "ThePEG/Persistency/PersistentOStream.h"
void MScale::persistentOutput(PersistentOStream & os) const {
  os <<theScaleChoice<<theMergingHelper;
}

#include "ThePEG/Persistency/PersistentIStream.h"
void MScale::persistentInput(PersistentIStream & is, int) {
  is >>theScaleChoice>>theMergingHelper;
}


  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<MScale,Herwig::MatchboxScaleChoice>
describeHerwigMScale("Herwig::MScale", "HwDipoleShower.so");

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
void MScale::Init() {
  
  static ClassDocumentation<MScale> documentation
  ("MScale implements scale choices related to transverse momenta.");
  
  static Reference<MScale,MatchboxScaleChoice> interfaceScaleChoice
  ("ScaleChoice",
   "Set the ScaleChoice to be considered.",
   &MScale::theScaleChoice, false, false, true, true, false);
  
  static Reference<MScale,Merger> interfaceMergingHelper
  ("MergingHelper",
   "",
   &MScale::theMergingHelper, false, false, true, true, false);
  
}

