// -*- C++ -*-
//
// MatchboxPhasespace.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxPhasespace class.
//

#include "MatchboxPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

MatchboxPhasespace::MatchboxPhasespace() {}

MatchboxPhasespace::~MatchboxPhasespace() {}

void MatchboxPhasespace::cloneDependencies(const std::string&) {
}

void MatchboxPhasespace::dumpInfo(const string& prefix) const {
  generator()->log() << prefix << fullName()
		     << " [" << this << "]\n";
  generator()->log() << prefix << "  | XComb " << lastXCombPtr()
		     << " for ";
  if ( lastXCombPtr() ) {
    for ( cPDVector::const_iterator p = lastXComb().mePartonData().begin();
	  p != lastXComb().mePartonData().end(); ++p ) {
      generator()->log() << (**p).PDGName() << " ";
    }
  }
  generator()->log() << "\n";
}

pair<double,Lorentz5Momentum> 
MatchboxPhasespace::timeLikeWeight(const Tree2toNDiagram& diag,
				   int branch) const {

  pair<int,int> children = diag.children(branch);

  if ( children.first == -1 ) {
    return make_pair(1.,meMomenta()[diag.externalId(branch)]);
  }

  pair<double,Lorentz5Momentum> res
    = timeLikeWeight(diag,children.first);

  pair<double,Lorentz5Momentum> other
    = timeLikeWeight(diag,children.second);

  res.first *= other.first;
  res.second += other.second;

  Energy2 mass2 = sqr(diag.allPartons()[branch]->mass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->width());

  res.first /=
    sqr((res.second.m2()-mass2)/lastSHat()) +
    mass2*width2/sqr(lastSHat());

  return res;

}

double MatchboxPhasespace::spaceLikeWeight(const Tree2toNDiagram& diag,
					   const Lorentz5Momentum& incoming,
					   int branch) const {

  if ( branch == -1 )
    return 1.;

  pair<int,int> children = diag.children(branch);

  pair<double,Lorentz5Momentum> res =
    timeLikeWeight(diag,children.second);

  if ( children.first == diag.nSpace() - 1 ) {
    return res.first;
  }

  res.second = incoming - res.second;

  Energy2 mass2 = sqr(diag.allPartons()[branch]->mass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->width());

  res.first /=
    sqr((res.second.m2()-mass2)/lastSHat()) +
    mass2*width2/sqr(lastSHat());

  return
    res.first * spaceLikeWeight(diag,res.second,children.first);

}

void MatchboxPhasespace::fillDiagramWeights() {

  diagramWeights.clear();

  for ( StandardXComb::DiagramVector::const_iterator d =
	  lastXComb().diagrams().begin(); d != lastXComb().diagrams().end(); ++d ) {
    diagramWeights[(**d).id()] = 
      spaceLikeWeight(dynamic_cast<const Tree2toNDiagram&>(**d),meMomenta()[0],0);
  }

}

Selector<MEBase::DiagramIndex> 
MatchboxPhasespace::selectDiagrams(const MEBase::DiagramVector& diags) const {
  Selector<MEBase::DiagramIndex> ret;
  for ( MEBase::DiagramIndex d = 0; d < diags.size(); ++d ) {
    ret.insert(diagramWeight(dynamic_cast<const Tree2toNDiagram&>(*diags[d])),d);
  }
  return ret;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).



// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractNoPIOClass<MatchboxPhasespace,HandlerBase>
  describeMatchboxPhasespace("Herwig::MatchboxPhasespace", "HwMatchbox.so");

void MatchboxPhasespace::Init() {

  static ClassDocumentation<MatchboxPhasespace> documentation
    ("MatchboxPhasespace defines an abstract interface to a phase "
     "space generator.");

}

