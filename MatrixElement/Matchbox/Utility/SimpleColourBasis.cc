// -*- C++ -*-
//
// SimpleColourBasis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

  if ( basis == id88 || basis == id33bar || basis == id33bar8 )
    return 1;

  if ( basis == id888 || basis == id33bar88 || basis == id33bar33bar )
    return 2;

  if ( basis == id8888 )
    return 6;

  throw Exception() << "SimpleColourBasis::prepareBasis(): Cannot handle colour configuration" << Exception::runerror;

  return 0;

}

double SimpleColourBasis::scalarProduct(size_t a, size_t b,
					const vector<PDT::Colour>& abBasis) const {

  if ( id33bar.empty() )
    makeIds();

  double Nc = SM().Nc();
  double Nc2 = sqr(Nc);
  double Nc3 = Nc*Nc2;
  double Nc4 = sqr(Nc2);
  double Nc6 = Nc2*Nc4;

  if ( a > b )
    swap(a,b);

  if ( !largeN() ) {

    if ( abBasis == id88 ) {
      return ( Nc2 - 1. )/4.;
    }

    if ( abBasis == id33bar ) {
      return Nc;
    }

    if ( abBasis == id888 ) {
      if ( a == b )
	return ( Nc4 - 3.*Nc2 + 2. )/(8.*Nc);
      return -( Nc2 - 1. )/(4.*Nc);
    }

    if ( abBasis == id33bar8 ) {
      return ( Nc2 - 1. )/2.;
    }

    if ( abBasis == id8888 ) {
      if ( a == b )
	return ( Nc6 - 4.*Nc4 + 6.*Nc2 - 3. )/(16.*Nc2);
      if ( ( a == 0 && b == 1 ) ||
	   ( a == 2 && b == 3 ) ||
	   ( a == 4 && b == 5 ) )
	return ( Nc4 + 2.*Nc2 - 3. )/(16.*Nc2);
      return -( Nc2 - 4. + 3./Nc2 )/16.;
    }

    if ( abBasis == id33bar88 ) {
      if ( a == b )
	return ( Nc4 - 2.*Nc2 + 1 )/(4.*Nc);
      return -( Nc2 - 1. )/(4.*Nc);
    }

    if ( abBasis == id33bar33bar ) {
      if ( a == b )
	return Nc2;
      return Nc;
    }

  } else {

    if ( a != b )
      return 0.;

    if ( abBasis == id88 ) {
      return Nc2/4.;
    }

    if ( abBasis == id33bar ) {
      return Nc;
    }

    if ( abBasis == id888 ) {
      return Nc3/8.;
    }

    if ( abBasis == id33bar8 ) {
      return Nc2/2.;
    }

    if ( abBasis == id8888 ) {
      return Nc4/16.;
    }

    if ( abBasis == id33bar88 ) {
      return Nc3/4.;
    }

    if ( abBasis == id33bar33bar ) {
      return Nc2;
    }

  }

  throw Exception() << "SimpleColourBasis::scalarProduct(): Cannot handle colour configuration" << Exception::runerror;

}

double SimpleColourBasis::tMatrixElement(size_t i, size_t a, 
					 size_t b,
					 const vector<PDT::Colour>&,
					 const vector<PDT::Colour>& bBasis) const {

  if ( id33bar.empty() )
    makeIds();

  if ( bBasis == id88 ) {
    if ( i == 0 )
      return a == 0 ? -1. : 1.;
    else
      return a == 0 ? 1. : -1.;
  }

  if ( bBasis == id33bar ) {
    return i == 0 ? 1. : -1.;
  }

  if ( bBasis == id888 ) {
    if ( i == 0 ) {
      if ( a == 3 && b == 0 )
	return 1.;
      if ( a == 0 && b == 0 )
	return -1.;
      if ( a == 1 && b == 1 )
	return 1.;
      if ( a == 2 && b == 1 )
	return -1.;
    }
    if ( i == 1 ) {
      if ( a == 4 && b == 0 )
	return 1.;
      if ( a == 3 && b == 0 )
	return -1.;
      if ( a == 2 && b == 1 )
	return 1.;
      if ( a == 5 && b == 1 )
	return -1.;
    }
    if ( i == 2 ) {
      if ( a == 0 && b == 0 )
	return 1.;
      if ( a == 4 && b == 0 )
	return -1.;
      if ( a == 5 && b == 1 )
	return 1.;
      if ( a == 1 && b == 1 )
	return -1.;
    }
    return 0.;
  }

  if ( bBasis == id33bar8 ) {
    if ( i == 0 )
      return a == 1 ? 1. : 0.;
    if ( i == 1 )
      return a == 0 ? -1. : 0.;
    if ( i == 2 )
      return a == 0 ? 1. : -1.;
  }

  throw Exception() << "SimpleColourBasis::tMatrixElement(): Cannot handle colour configuration" << Exception::runerror;

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

  if ( basis == id88 ) {
    return
      ( idColoured == 0 && idAntiColoured == 1 ) ||
      ( idColoured == 1 && idAntiColoured == 0 );
  }

  if ( basis == id33bar ) {
    return
      idColoured == 0 && idAntiColoured == 1;
  }

  if ( basis == id888 ) {
    if ( a == 0 )
      return
	( idColoured == 0 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 0 );
    if ( a == 1 )
      return
	( idColoured == 0 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 0 );
  }

  if ( basis == id33bar8 ) {
    return
      ( idColoured == 0 && idAntiColoured == 2 ) ||
      ( idColoured == 2 && idAntiColoured == 1 );
  }

  if ( basis == id8888 ) {
    if ( a == 0 )
      return
	( idColoured == 0 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 0 );
    if ( a == 1 )
      return
	( idColoured == 0 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 0 );
    if ( a == 2 )
      return
	( idColoured == 0 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 0 );
    if ( a == 3 )
      return
	( idColoured == 0 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 0 );
    if ( a == 4 )
      return
	( idColoured == 0 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 0 );
    if ( a == 5 )
      return
	( idColoured == 0 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 1 ) ||
	( idColoured == 1 && idAntiColoured == 0 );
  }

  if ( basis == id33bar88 ) {
    if ( a == 0 )
      return
	( idColoured == 0 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 1 );
    if ( a == 1 )
      return
	( idColoured == 0 && idAntiColoured == 3 ) ||
	( idColoured == 3 && idAntiColoured == 2 ) ||
	( idColoured == 2 && idAntiColoured == 1 );
  }

  if ( basis == id33bar33bar ) {
    if ( a == 0 )
      return
	( idColoured == 0 && idAntiColoured == 1 ) ||
	( idColoured == 2 && idAntiColoured == 3 );
    if ( a == 1 )
      return
	( idColoured == 0 && idAntiColoured == 3 ) ||
	( idColoured == 2 && idAntiColoured == 1 );
  }

  return false;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void SimpleColourBasis::makeIds() const {

  id88.push_back(PDT::Colour8);
  id88.push_back(PDT::Colour8);

  id33bar.push_back(PDT::Colour3);
  id33bar.push_back(PDT::Colour3bar);

  id888.push_back(PDT::Colour8);
  id888.push_back(PDT::Colour8);
  id888.push_back(PDT::Colour8);

  id33bar8.push_back(PDT::Colour3);
  id33bar8.push_back(PDT::Colour3bar);
  id33bar8.push_back(PDT::Colour8);

  id8888.push_back(PDT::Colour8);
  id8888.push_back(PDT::Colour8);
  id8888.push_back(PDT::Colour8);
  id8888.push_back(PDT::Colour8);

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
  describeHerwigSimpleColourBasis("Herwig::SimpleColourBasis", "Herwig.so");

void SimpleColourBasis::Init() {

  static ClassDocumentation<SimpleColourBasis> documentation
    ("SimpleColourBasis implements the colour algebra needed for "
     "electroweak boson and electroweak boson + jet production at NLO. It mainly "
     "serves as an example for the general ColourBasis interface.");

}

