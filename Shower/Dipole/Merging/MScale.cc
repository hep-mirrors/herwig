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

#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

#include "ThePEG/Cuts/JetFinder.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MScale::MScale() {}

MScale::~MScale() {}

IBPtr MScale::clone() const {
  return new_ptr(*this);
}

IBPtr MScale::fullclone() const {
  return new_ptr(*this);
}


Energy2 MScale::renormalizationScale() const {

  
  if (theMergingHelper->renormscale()/GeV> Constants::epsilon ){
      return sqr(theMergingHelper->renormscale());
  }
  if (!theScaleChoice){
    cerr<<"MScale: You need to set a scale for MScale.";
  }
  theScaleChoice->setXComb(theLastXComb);
  return theScaleChoice->renormalizationScale();
}

Energy2 MScale::factorizationScale() const {
  if (theMergingHelper->renormscale()/GeV> Constants::epsilon )
      return sqr( theMergingHelper->renormscale());
  
  return sqr(10.*GeV);
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MScale::persistentOutput(PersistentOStream & os) const {
  os <<theScaleChoice<<theMergingHelper;
}

void MScale::persistentInput(PersistentIStream & is, int) {
  is >>theScaleChoice>>theMergingHelper;
}


  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MScale,Herwig::MatchboxScaleChoice>
describeHerwigMScale("Herwig::MScale", "HwDipoleShower.so");

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

