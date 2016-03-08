  // -*- C++ -*-
  //
  // MatchboxMergeScale.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2012 The Herwig Collaboration
  //
  // Herwig is licenced under version 2 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
  //
  // This is the implementation of the non-inlined, non-templated member
  // functions of the MatchboxMergeScale class.
  //

#include "MatchboxMergeScale.h"
#include "ClusterNode.h"

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

MatchboxMergeScale::MatchboxMergeScale() {}

MatchboxMergeScale::~MatchboxMergeScale() {}

IBPtr MatchboxMergeScale::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMergeScale::fullclone() const {
  return new_ptr(*this);
}


Energy2 MatchboxMergeScale::renormalizationScale() const {
    //cout<<"\nrenscale";
  Ptr<MatchboxMEBase>::ptr ME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(lastXComb().matrixElement());
  if (ME)
    if (ME->firstNode()->renormscale()/GeV> Constants::epsilon ){
        //cout<<"\n"<<ME->firstNode()->renormscale()/GeV;
      return sqr(ME->firstNode()->renormscale());
    }
  if (!theScaleChoice){
    cerr<<"MatchboxMergeScale: You need to set a scale for MatchboxMergeScale.";
  }
  theScaleChoice->setXComb(theLastXComb);
    //cout<<"\n-->"<<theScaleChoice->renormalizationScale()/GeV2;
  return theScaleChoice->renormalizationScale();
  
  
  
}

Energy2 MatchboxMergeScale::factorizationScale() const {
  Ptr<MatchboxMEBase>::ptr ME = dynamic_ptr_cast<Ptr<MatchboxMEBase>::ptr>(lastXComb().matrixElement());
  if (ME)
    if (ME->firstNode()->renormscale()/GeV> Constants::epsilon )
      return sqr( ME->firstNode()->renormscale());
  
  return sqr(10.*GeV);
}

  // If needed, insert default implementations of virtual function defined
  // in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxMergeScale::persistentOutput(PersistentOStream & os) const {
  os <<theScaleChoice;
}

void MatchboxMergeScale::persistentInput(PersistentIStream & is, int) {
  is >>theScaleChoice;
}


  // *** Attention *** The following static variable is needed for the type
  // description system in ThePEG. Please check that the template arguments
  // are correct (the class and its base class), and that the constructor
  // arguments are correct (the class name and the name of the dynamically
  // loadable library where the class implementation can be found).
DescribeClass<MatchboxMergeScale,MatchboxScaleChoice>
describeHerwigMatchboxMergeScale("Herwig::MatchboxMergeScale", "Herwig.so");

void MatchboxMergeScale::Init() {
  
  static ClassDocumentation<MatchboxMergeScale> documentation
  ("MatchboxMergeScale implements scale choices related to transverse momenta.");
  
  
  static Reference<MatchboxMergeScale,MatchboxScaleChoice> interfaceScaleChoice
  ("ScaleChoice",
   "Set the ScaleChoice to be considered.",
   &MatchboxMergeScale::theScaleChoice, false, false, true, true, false);
  
  
  
  
}

