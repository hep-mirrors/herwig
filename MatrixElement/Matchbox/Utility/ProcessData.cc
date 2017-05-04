// -*- C++ -*-
//
//ProcessData.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProcessData class.
//

#include "ProcessData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ProcessData::ProcessData() {}

ProcessData::~ProcessData() {
  for ( map<Ptr<Tree2toNDiagram>::tcptr,vector<ColourLines*> >::iterator cl =
	  theColourFlowMap.begin(); cl != theColourFlowMap.end(); ++cl ) {
    for ( vector<ColourLines*>::iterator c = cl->second.begin();
	  c != cl->second.end(); ++c ) {
      if ( *c )
	delete *c;
    }
  }
  theColourFlowMap.clear();
}

IBPtr ProcessData::clone() const {
  return new_ptr(*this);
}

IBPtr ProcessData::fullclone() const {
  return new_ptr(*this);
}

void ProcessData::fillMassGenerators(const PDVector& subpro) {
  for ( PDVector::const_iterator pd = subpro.begin();
	pd != subpro.end(); ++pd ) {
    if ( theMassGenerators.find(*pd) != theMassGenerators.end() )
      continue;
    tGenericMassGeneratorPtr mgen = 
      dynamic_ptr_cast<tGenericMassGeneratorPtr>((**pd).massGenerator());
    if ( !mgen )
      continue;
    theMassGenerators[*pd] = mgen;
  }
}

tGenericMassGeneratorPtr ProcessData::massGenerator(cPDPtr pd) {
  if ( !pd->massGenerator() )
    return tGenericMassGeneratorPtr();
  map<cPDPtr,tGenericMassGeneratorPtr>::iterator mg =
    theMassGenerators.find(pd);
  if ( mg == theMassGenerators.end() )
    return tGenericMassGeneratorPtr();
  return mg->second;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ProcessData::persistentOutput(PersistentOStream & os) const {
  os << theDiagramMap << theMassGenerators;
}

void ProcessData::persistentInput(PersistentIStream & is, int) {
  is >> theDiagramMap >> theMassGenerators;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<ProcessData,HandlerBase>
  describeHerwigProcessData("Herwig::ProcessData", "Herwig.so");

void ProcessData::Init() {

  static ClassDocumentation<ProcessData> documentation
    ("Provide storage for process data");

}

