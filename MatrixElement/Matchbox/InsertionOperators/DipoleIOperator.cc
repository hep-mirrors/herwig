// -*- C++ -*-
//
// DipoleIOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleIOperator class.
//

#include "DipoleIOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"

using namespace Herwig;
using Constants::pi;

DipoleIOperator::DipoleIOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    betaZero(-1.),
    KQuark(-1.0), KGluon(-1.0) {}

DipoleIOperator::~DipoleIOperator() {}

IBPtr DipoleIOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleIOperator::fullclone() const {
  return new_ptr(*this);
}

void DipoleIOperator::setBorn(Ptr<MatchboxMEBase>::tptr me) {
  MatchboxInsertionOperator::setBorn(me);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLight();
    betaZero = gammaGluon;
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLight();
    if ( isDR() ) {
      gammaQuark -= CF/2.;
      gammaGluon -= CA/6.;
    }
  }
}

bool DipoleIOperator::apply(const cPDVector& pd) const {
  bool first = false;
  bool second = false;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {
    if ( !first ) {
      if ( apply(*p) )
	first = true;
    } else {
      if ( apply(*p) )
	second = true;
    }
  }
  return first && second;
}

bool DipoleIOperator::apply(tcPDPtr pd) const {
  return
    pd->mass() == ZERO &&
    (abs(pd->id()) < 6 || pd->id() == ParticleID::g);
}

double DipoleIOperator::me2() const {

  double res = 0.;

  int idi = 0; int idj = 0;

  Energy2 mu2 = lastBorn()->mu2();

  for ( cPDVector::const_iterator i = mePartonData().begin();
	i != mePartonData().end(); ++i, ++idi ) {

    if ( !apply(*i) )
      continue;

    idj = 0;

    for ( cPDVector::const_iterator j = mePartonData().begin();
	  j != mePartonData().end(); ++j, ++idj ) {

      if ( !apply(*j) )
	continue;

      if ( i == j || lastBorn()->noDipole(idi,idj) )
	continue;

      double delta = 0.;

      double xgammaGluon = gammaGluon;
      double xgammaQuark = gammaQuark;
      if ( isDR() ) {
	xgammaGluon += CA/6.;
	xgammaQuark += CF/2.;
      }

      if ( isBDK() ) {
	assert(!isCS() && !isExpanded());
	delta = 
	  ((**i).id() == ParticleID::g ? xgammaGluon : xgammaQuark) * 
	  log(mu2/(2.*(meMomenta()[idi]*meMomenta()[idj])));
	if ( (idi < 2 && idj < 2) ||
	     (idi > 1 && idj > 1) )
	  delta +=
	    ((**i).id() == ParticleID::g ? CA : CF) * sqr(pi) / 2.;
      }

      if ( isExpanded() ) {
	assert(!isCS() && !isBDK());
	double theLog = log(mu2/(2.*(meMomenta()[idi]*meMomenta()[idj])));
	delta = 
	  ((**i).id() == ParticleID::g ? CA : CF) * 0.5 * sqr(theLog) +
	  ((**i).id() == ParticleID::g ? xgammaGluon : xgammaQuark) * theLog;
      }

      res +=
	( ((**i).id() == ParticleID::g ? CA : CF) * (-sqr(pi)/3.) +
	  ((**i).id() == ParticleID::g ? gammaGluon : gammaQuark) +
	  ((**i).id() == ParticleID::g ? KGluon : KQuark) +
	  delta ) *
	lastBorn()->colourCorrelatedME2(make_pair(idi,idj));
    }
  }

  Energy2 muR2 = 
    lastBorn()->renormalizationScale()*
    lastBorn()->renormalizationScaleFactor();
  if ( muR2 != mu2 ) {
    res -=
      betaZero *
      lastBorn()->orderInAlphaS() * log(muR2/mu2) *
      lastBorn()->me2();
  }

  // include the finite renormalization for DR here; ATTENTION this
  // has to be mentioned in the manual!  see hep-ph/9305239 for
  // details; this guarantees an expansion in alpha_s^\bar{MS} when
  // using dimensional reduction
  if ( isDR() )
    res -= (CA/6.)*lastBorn()->orderInAlphaS()*lastBorn()->me2();

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );

  return res;

}

double DipoleIOperator::oneLoopDoublePole() const {

  if ( !isExpanded() )
    return 0.;

  double res = 0.;

  for ( cPDVector::const_iterator i = mePartonData().begin();
	i != mePartonData().end(); ++i ) {

    if ( !apply(*i) )
      continue;
    
    res += (**i).id() == ParticleID::g ? CA : CF;

  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) ) * ( - lastBorn()->me2() );

  return res;

}

double DipoleIOperator::oneLoopSinglePole() const {

  if ( !isExpanded() )
    return 0.;

  double res = 0.;

  int idi = 0; int idj = 0;

  Energy2 mu2 = lastBorn()->mu2();

  for ( cPDVector::const_iterator i = mePartonData().begin();
	i != mePartonData().end(); ++i, ++idi ) {

    if ( !apply(*i) )
      continue;

    idj = 0;

    for ( cPDVector::const_iterator j = mePartonData().begin();
	  j != mePartonData().end(); ++j, ++idj ) {

      if ( !apply(*j) )
	continue;

      if ( i == j || lastBorn()->noDipole(idi,idj) )
	continue;

      double xgammaGluon = gammaGluon;
      double xgammaQuark = gammaQuark;
      if ( isDR() ) {
	xgammaGluon += CA/6.;
	xgammaQuark += CF/2.;
      }

      double theLog = log(mu2/(2.*(meMomenta()[idi]*meMomenta()[idj])));
      double delta = 
	((**i).id() == ParticleID::g ? CA : CF) * theLog +
	((**i).id() == ParticleID::g ? xgammaGluon : xgammaQuark);

      res +=
	delta * lastBorn()->colourCorrelatedME2(make_pair(idi,idj));
    }

  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );

  return res;

}

void DipoleIOperator::persistentOutput(PersistentOStream & os) const {
  os << CA << CF << gammaQuark << gammaGluon << betaZero
     << KQuark << KGluon;
}

void DipoleIOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> betaZero
     >> KQuark >> KGluon;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleIOperator,MatchboxInsertionOperator>
describeHerwigDipoleIOperator("Herwig::DipoleIOperator", "HwMatchbox.so");

void DipoleIOperator::Init() {

  static ClassDocumentation<DipoleIOperator> documentation
    ("DipoleIOperator");

  DipoleRepository::registerInsertionOperator<DipoleIOperator>("LightIOperator");

}

