// -*- C++ -*-
//
// MatchboxHybridAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxHybridAmplitude class.
//

#include "MatchboxHybridAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/MatchboxMEBase.h"

using namespace Herwig;

MatchboxHybridAmplitude::MatchboxHybridAmplitude() 
  : theUseOLPCorrelators(false) {}

MatchboxHybridAmplitude::~MatchboxHybridAmplitude() {}

IBPtr MatchboxHybridAmplitude::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxHybridAmplitude::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxHybridAmplitude::isConsistent() const {
  assert(oneLoopAmplitude());
  return
    !treeLevelAmplitude()->isOLPTree() &&
    !treeLevelAmplitude()->isOLPLoop() &&
    oneLoopAmplitude()->haveOneLoop() &&
    treeLevelAmplitude()->orderInGs() ==
    oneLoopAmplitude()->orderInGs() &&
    treeLevelAmplitude()->orderInGem() ==
    oneLoopAmplitude()->orderInGem() &&
    treeLevelAmplitude()->hasRunningAlphaS() ==
    oneLoopAmplitude()->hasRunningAlphaS() &&
    treeLevelAmplitude()->hasRunningAlphaEW() ==
    oneLoopAmplitude()->hasRunningAlphaEW() &&
    !(treeLevelAmplitude()->nDimAdditional() != 0 &&
      oneLoopAmplitude()->nDimAdditional() != 0);
}

bool MatchboxHybridAmplitude::canHandle(const PDVector& p,
					Ptr<MatchboxFactory>::tptr f,
					bool virt) const { 
  if ( !virt )
    return treeLevelAmplitude()->canHandle(p,f,false);
  if ( treeLevelAmplitude()->canHandle(p,f,false) &&
       oneLoopAmplitude()->canHandle(p,f,true) ) {
    if ( !isConsistent() ) {
      generator()->log() << "Warning: Inconsistent settings encountered for MatchboxHybridAmplitude '"
			 << name() << "'\n" << flush;
      return false;
    }
    return true;
  }
  return false;
}

void MatchboxHybridAmplitude::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  treeLevelAmplitude()->prepareAmplitudes(me);
}

void MatchboxHybridAmplitude::prepareOneLoopAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  assert(oneLoopAmplitude());
  oneLoopAmplitude()->prepareOneLoopAmplitudes(me);
}

double MatchboxHybridAmplitude::symmetryRatio() const {

  assert(oneLoopAmplitude());

  double ifact = 1.;
  if ( treeLevelAmplitude()->hasInitialAverage() &&
       !oneLoopAmplitude()->hasInitialAverage() ) {
    ifact = 1./4.;
    if (lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour3 ||
	 lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour3bar )
      ifact /= SM().Nc();
    else if ( lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour8 )
      ifact /= (SM().Nc()*SM().Nc()-1.);

    if ( lastMatchboxXComb()->matchboxME()->mePartonData()[1]->iColour() == PDT::Colour3 ||
	 lastMatchboxXComb()->matchboxME()->mePartonData()[1]->iColour() == PDT::Colour3bar )
      ifact /= SM().Nc();
    else if ( mePartonData()[1]->iColour() == PDT::Colour8 )
      ifact /= (SM().Nc()*SM().Nc()-1.);
  }

  if ( !treeLevelAmplitude()->hasInitialAverage() &&
       oneLoopAmplitude()->hasInitialAverage() ) {
    ifact = 4.;
    if ( lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour3 ||
	 lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour3bar )
      ifact *= SM().Nc();
    else if ( lastMatchboxXComb()->matchboxME()->mePartonData()[0]->iColour() == PDT::Colour8 )
      ifact *= (SM().Nc()*SM().Nc()-1.);

    if ( lastMatchboxXComb()->matchboxME()->mePartonData()[1]->iColour() == PDT::Colour3 ||
	 lastMatchboxXComb()->matchboxME()->mePartonData()[1]->iColour() == PDT::Colour3bar )
      ifact *= SM().Nc();
    else if ( lastMatchboxXComb()->matchboxME()->mePartonData()[1]->iColour() == PDT::Colour8 )
      ifact *= (SM().Nc()*SM().Nc()-1.);
  }

  if ( treeLevelAmplitude()->hasFinalStateSymmetry() &&
       !oneLoopAmplitude()->hasFinalStateSymmetry() ) {
    assert(lastMatchboxXComb()->matchboxME());
    ifact *= lastMatchboxXComb()->matchboxME()->finalStateSymmetry();
  }

  if ( !treeLevelAmplitude()->hasFinalStateSymmetry() &&
       oneLoopAmplitude()->hasFinalStateSymmetry() ) {
    assert(lastMatchboxXComb()->matchboxME());
    ifact /= lastMatchboxXComb()->matchboxME()->finalStateSymmetry();
  }

  return ifact;

}

void MatchboxHybridAmplitude::cloneDependencies(const std::string& prefix,bool slim) {

  if ( treeLevelAmplitude() ) {
    Ptr<MatchboxAmplitude>::ptr myTreeLevelAmplitude = treeLevelAmplitude()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myTreeLevelAmplitude->name();
    if ( ! (generator()->preinitRegister(myTreeLevelAmplitude,pname.str()) ) )
      throw Exception() << "MatchboxHybridAmplitude::cloneDependencies(): Amplitude " << pname.str() << " already existing." << Exception::runerror;
    myTreeLevelAmplitude->cloneDependencies(pname.str());
    treeLevelAmplitude(myTreeLevelAmplitude);
  }

  if ( oneLoopAmplitude() ) {
    Ptr<MatchboxAmplitude>::ptr myOneLoopAmplitude = oneLoopAmplitude()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myOneLoopAmplitude->name();
    if ( ! (generator()->preinitRegister(myOneLoopAmplitude,pname.str()) ) )
      throw Exception() << "MatchboxHybridAmplitude::cloneDependencies(): Amplitude " << pname.str() << " already existing." << Exception::runerror;
    myOneLoopAmplitude->cloneDependencies(pname.str());
    oneLoopAmplitude(myOneLoopAmplitude);
  }

  MatchboxAmplitude::cloneDependencies(prefix,slim);

}

void MatchboxHybridAmplitude::doinit() {
  MatchboxAmplitude::doinit();
  if ( treeLevelAmplitude() )
    treeLevelAmplitude()->init();
  if ( oneLoopAmplitude() )
    oneLoopAmplitude()->init();
}

void MatchboxHybridAmplitude::doinitrun() {
  MatchboxAmplitude::doinitrun();
  if ( treeLevelAmplitude() )
    treeLevelAmplitude()->initrun();
  if ( oneLoopAmplitude() )
    oneLoopAmplitude()->initrun();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxHybridAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theTreeLevelAmplitude << theOneLoopAmplitude << theUseOLPCorrelators;
}

void MatchboxHybridAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theTreeLevelAmplitude >> theOneLoopAmplitude >> theUseOLPCorrelators;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxHybridAmplitude,Herwig::MatchboxAmplitude>
  describeHerwigMatchboxHybridAmplitude("Herwig::MatchboxHybridAmplitude", "Herwig.so");

void MatchboxHybridAmplitude::Init() {

  static ClassDocumentation<MatchboxHybridAmplitude> documentation
    ("MatchboxHybridAmplitude unifies two amplitude objects to "
     "provide tree and one-loop matrix elements.");


  static Reference<MatchboxHybridAmplitude,MatchboxAmplitude> interfaceTreeLevelAmplitude
    ("TreeLevelAmplitude",
     "Set the tree level amplitude to be used.",
     &MatchboxHybridAmplitude::theTreeLevelAmplitude, false, false, true, false, false);

  static Reference<MatchboxHybridAmplitude,MatchboxAmplitude> interfaceOneLoopAmplitude
    ("OneLoopAmplitude",
     "Set the one-loop amplitude to be used.",
     &MatchboxHybridAmplitude::theOneLoopAmplitude, false, false, true, true, false);

  static Switch<MatchboxHybridAmplitude,bool> interfaceUseOLPCorrelators
    ("UseOLPCorrelators",
     "Obtain correlated matrix elements from the OLP instead of "
     "the tree-level amplitude.",
     &MatchboxHybridAmplitude::theUseOLPCorrelators, false, false, false);
  static SwitchOption interfaceUseOLPCorrelatorsYes
    (interfaceUseOLPCorrelators,
     "Yes",
     "",
     true);
  static SwitchOption interfaceUseOLPCorrelatorsNo
    (interfaceUseOLPCorrelators,
     "No",
     "",
     false);
  interfaceUseOLPCorrelators.rank(-1);

}

