// -*- C++ -*-
//
// MatchboxReference.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxReference class.
//

#include "MatchboxReference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxReference::MatchboxReference() {}

MatchboxReference::~MatchboxReference() {}

IBPtr MatchboxReference::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxReference::fullclone() const {
  return new_ptr(*this);
}

double MatchboxReference::generateTwoToNKinematics(const double*,
						   vector<Lorentz5Momentum>& momenta) {


  map<cPDVector,ifstream*>::iterator ref =
    referenceSamples.find(mePartonData());
  if ( ref == referenceSamples.end() ) {
    ostringstream refname;
    for ( cPDVector::const_iterator p = mePartonData().begin();
	  p != mePartonData().end(); ++p ) {
      refname << (**p).PDGName();
    }
    refname << ".rambo";
    referenceSamples[mePartonData()] = new ifstream(refname.str().c_str());
    ref = referenceSamples.find(mePartonData());
  }
  assert(ref != referenceSamples.end());

  ifstream& in = *(ref->second);
  assert(in);

  double x1,x2;
  double x,y,z,t,m;
  double weight;

  in >> x1 >> x2;
  for ( vector<Lorentz5Momentum>::iterator p = momenta.begin();
	p != momenta.end(); ++p ) {
    in >> x >> y >> z >> t >> m;
    *p = Lorentz5Momentum(x*GeV,y*GeV,z*GeV,t*GeV,m*GeV);
  }

  in >> weight;

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  fillDiagramWeights();

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxReference::persistentOutput(PersistentOStream &) const {}

void MatchboxReference::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxReference,MatchboxPhasespace>
  describeHerwigMatchboxReference("Herwig::MatchboxReference", "Herwig.so");

void MatchboxReference::Init() {

  static ClassDocumentation<MatchboxReference> documentation
    ("MatchboxReference implements reference sample phase space generation.");

}

