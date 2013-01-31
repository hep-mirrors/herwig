// -*- C++ -*-
//
// PowhegRealReweight.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PowhegRealReweight class.
//

#include "PowhegRealReweight.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PowhegRealReweight::PowhegRealReweight()
  : ME2byDipoles(), theBornScreening(true) {}

PowhegRealReweight::~PowhegRealReweight() {}

IBPtr PowhegRealReweight::clone() const {
  return new_ptr(*this);
}

IBPtr PowhegRealReweight::fullclone() const {
  return new_ptr(*this);
}

double PowhegRealReweight::evaluate() const {

  if ( !bornScreening() )
    return 0.0;

  double ratio = ME2byDipoles::evaluate();

  double born = scaledBorn();
  double screen = scaledBornScreen();

  ratio *= screen / ( born + screen );

  return ratio;

}

void PowhegRealReweight::setVetoScales(tSubProPtr subpro) const {

  Energy pt =
    projectionDipole()->lastPt();

  if ( projectionDipole()->realEmitter() == 0 ||
       projectionDipole()->realSpectator() == 0 ) {
    if ( subpro->incoming().first->vetoScale() < 0.0*GeV2 ||
	 subpro->incoming().first->vetoScale() > sqr(pt) )
      subpro->incoming().first->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realEmitter() == 1 ||
       projectionDipole()->realSpectator() == 1 ) {
    if ( subpro->incoming().second->vetoScale() < 0.0*GeV2 ||
	 subpro->incoming().second->vetoScale() > sqr(pt) )
      subpro->incoming().second->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realEmitter() > 1 ) {
    if ( subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale() < 0.0*GeV2 ||
	 subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale() > sqr(pt) )
      subpro->outgoing()[projectionDipole()->realEmitter()-2]->vetoScale(sqr(pt));
  }

  if ( projectionDipole()->realSpectator() > 1 ) {
    if ( subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale() < 0.0*GeV2 ||
	 subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale() > sqr(pt) )
      subpro->outgoing()[projectionDipole()->realSpectator()-2]->vetoScale(sqr(pt));
  }

  if ( subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale() < 0.0*GeV2 ||
       subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale() > sqr(pt) )
    subpro->outgoing()[projectionDipole()->realEmission()-2]->vetoScale(sqr(pt));  

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PowhegRealReweight::persistentOutput(PersistentOStream & os) const {
  os << theBornScreening;
}

void PowhegRealReweight::persistentInput(PersistentIStream & is, int) {
  is >> theBornScreening;
}

void PowhegRealReweight::Init() {

  static ClassDocumentation<PowhegRealReweight> documentation
    ("PowhegRealReweight");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PowhegRealReweight,ME2byDipoles>
describeHerwigPowhegRealReweight("Herwig::PowhegRealReweight", "HwMatchbox.so");
