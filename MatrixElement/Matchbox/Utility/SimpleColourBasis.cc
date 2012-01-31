// -*- C++ -*-
//
// SimpleColourBasis.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleColourBasis class.
//

#include "SimpleColourBasis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SimpleColourBasis::SimpleColourBasis() {}

SimpleColourBasis::~SimpleColourBasis() {}

IBPtr SimpleColourBasis::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleColourBasis::fullclone() const {
  return new_ptr(*this);
}

size_t SimpleColourBasis::prepareBasis(const vector<PDT::Colour>& basis) {

  if ( id33bar.empty() )
    makeIds();

  if ( basis == id33bar || basis == id33bar8 )
    return 1;

  if ( basis == id33bar33bar )
    return 2;

  if ( basis == id33bar88 )
    return 2;

  throw Exception() << "Cannot handle colour configuration" << Exception::abortnow;

  return 0;

}

double SimpleColourBasis::scalarProduct(size_t a, size_t b,
					const vector<PDT::Colour>& abBasis) const {

  if ( id33bar.empty() )
    makeIds();

  if ( abBasis == id33bar ) {
    assert(a==b);
    return SM().Nc();
  }

  if ( abBasis == id33bar8 ) {
    assert(a==b);
    return (SM().Nc()*SM().Nc()-1.)/2.;
  }

  if ( abBasis == id33bar88 ) {
    if ( a == b ) {
      return sqr((SM().Nc()*SM().Nc()-1.)/2.)/SM().Nc();
    }
    return 1./(4.*SM().Nc());
  }

  if ( abBasis == id33bar33bar ) {
    if ( a == b )
      return SM().Nc()*SM().Nc();
    return SM().Nc();
  }

  throw Exception() << "Cannot handle colour configuration" << Exception::abortnow;

}

double SimpleColourBasis::tMatrixElement(size_t i, size_t a, size_t b,
					 const vector<PDT::Colour>&,
					 const vector<PDT::Colour>& bBasis) const {

  if ( id33bar.empty() )
    makeIds();

  if ( bBasis == id33bar ) {
    assert(a==b);
    return i == 0 ? 1. : -1.;
  }

  if ( bBasis == id33bar8 ) {
    assert(b==0);
    if ( i == 0 )
      return a == 0 ? 1. : 0.;
    if ( i == 1 )
      return a == 1 ? -1. : 0.;
    if ( i == 2 && a != 2 )
      return a == 0 ? -1. : 1.;
    return 0.;
  }

  throw Exception() << "Cannot handle colour configuration" << Exception::abortnow;

  return 0.;

}

bool SimpleColourBasis::colourConnected(const cPDVector& sub,
					const vector<PDT::Colour>& basis,
					const pair<int,bool>& i, 
					const pair<int,bool>& j, 
					size_t a) const {

  if ( id33bar.empty() )
    makeIds();

  // translate process to basis ids
  map<cPDVector,map<size_t,size_t> >::const_iterator trans
    = indexMap().find(sub);
  assert(trans != indexMap().end());

  int idColoured = i.second ? j.first : i.first;
  idColoured = trans->second.find(idColoured)->second;
  int idAntiColoured = i.second ? i.first : j.first;
  idAntiColoured = trans->second.find(idAntiColoured)->second;

  if ( basis == id33bar )
    return idColoured == 0 && idAntiColoured == 1;

  if ( basis == id33bar8 ) {
    return 
      (idColoured == 0 && idAntiColoured == 2) || 
      (idColoured == 2 && idAntiColoured == 1);
  }

  if ( basis == id33bar88 ) {
    if ( a == 0 )
      return
	(idColoured == 0 && idAntiColoured == 2) || 
	(idColoured == 3 && idAntiColoured == 1) ||
	(idColoured == 2 && idAntiColoured == 3);
    if ( a == 1 )
      return
	(idColoured == 0 && idAntiColoured == 3) || 
	(idColoured == 2 && idAntiColoured == 1) ||
	(idColoured == 3 && idAntiColoured == 2);
  }

  if ( basis == id33bar33bar ) {
    if ( a == 0 )
      return
	(idColoured == 0 && idAntiColoured == 1) ||
	(idColoured == 2 && idAntiColoured == 3);
    if ( a == 1 )
      return
	(idColoured == 0 && idAntiColoured == 3) ||
	(idColoured == 2 && idAntiColoured == 1);
  }

  return false;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void SimpleColourBasis::makeIds() const {

  id33bar.push_back(PDT::Colour3);
  id33bar.push_back(PDT::Colour3bar);

  id33bar8.push_back(PDT::Colour3);
  id33bar8.push_back(PDT::Colour3bar);
  id33bar8.push_back(PDT::Colour8);

  id33bar88.push_back(PDT::Colour3);
  id33bar88.push_back(PDT::Colour3bar);
  id33bar88.push_back(PDT::Colour8);
  id33bar88.push_back(PDT::Colour8);

  id33bar33bar.push_back(PDT::Colour3);
  id33bar33bar.push_back(PDT::Colour3bar);
  id33bar33bar.push_back(PDT::Colour3);
  id33bar33bar.push_back(PDT::Colour3bar);

}

void SimpleColourBasis::persistentOutput(PersistentOStream &) const {}

void SimpleColourBasis::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<SimpleColourBasis,ColourBasis>
  describeHerwigSimpleColourBasis("Herwig::SimpleColourBasis", "HwMatchbox.so");

void SimpleColourBasis::Init() {

  static ClassDocumentation<SimpleColourBasis> documentation
    ("SimpleColourBasis implements the colour algebra needed for "
     "vector boson and vector boson + jet production at NLO. It mainly "
     "serves as an example for the general ColourBasis interface.");

}

