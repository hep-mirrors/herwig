// -*- C++ -*-
//
// DipoleIOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

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

//////////////////////////////////////////////////////////////////////

bool DipoleIOperator::apply(const cPDVector& pd) const {

  // DipoleIOperator should only apply if in the overall
  // process only massless partons can occur.

  // Prohibit splittings g->Q\bar{Q} in the final state.
  // These are covered completely by DipoleMIOperator.
  if ( NHeavyJetVec().size()!=0 ) {
    return false;
  }

  bool first = false;
  bool second = false;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {
    // Since this loop only checks for at least one exis-
    // ting combination: Return false if any massive par-
    // tons are present (covered by DipoleMIOperator).
    if ( (*p)->coloured() && (*p)->hardProcessMass()!=ZERO ) {
      return false;
    }
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
    pd->hardProcessMass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

void DipoleIOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
      // gammaGluon = (11./6.)*CA - (1./3.)*NLightJetVec().size();
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLightJetVec().size();
    betaZero = gammaGluon;
    double alpha = factory()->alphaParameter();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    KQuark +=-CF*sqr(log(alpha))+gammaQuark*(alpha-1-log(alpha));
    
    
      // KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*NLightJetVec().size();
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLightJetVec().size();
    KGluon +=-CA*sqr(log(alpha))+gammaGluon*(alpha-1-log(alpha));
                                            
                                            
    if ( isDR() ) {
      gammaQuark -= CF/2.;
      gammaGluon -= CA/6.;
    }
  }
}

//////////////////////////////////////////////////////////////////////

vector<int> DipoleIOperator::NLightJetVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleIOperator::NLightJetVec(): Could not find a jet particle group named 'j'" << Exception::runerror;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNLightJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).hardProcessMass() == ZERO )
      theNLightJetVec.push_back( (**theP).id() );
  }

  return theNLightJetVec;

}

vector<int> DipoleIOperator::NHeavyJetVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleIOperator::NHeavyJetVec(): Could not find a jet particle group named 'j'" << Exception::runerror;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNHeavyJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).hardProcessMass() != ZERO )
      theNHeavyJetVec.push_back( (**theP).id() );
  }

  return theNHeavyJetVec;

}

vector<int> DipoleIOperator::NLightBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNLightBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    if ( abs((**j).id()) < 7 && (**j).hardProcessMass() == ZERO )
      theNLightBornVec.push_back( (**j).id() );
  }

  return theNLightBornVec;

}

vector<int> DipoleIOperator::NHeavyBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNHeavyBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    if ( abs((**j).id()) < 7 && (**j).hardProcessMass() != ZERO )
      theNHeavyBornVec.push_back( (**j).id() );
  }

  return theNHeavyBornVec;

}

vector<int> DipoleIOperator::NLightProtonVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("p");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleIOperator::NLightProtonVec(): Could not find a proton particle group named 'p'" << Exception::runerror;

  const PDVector& theProtonConstitutents = theIt->second;
  vector<int> theNLightProtonVec;

  for ( PDVector::const_iterator theP = theProtonConstitutents.begin();
        theP != theProtonConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).hardProcessMass() == ZERO )
      theNLightProtonVec.push_back( (**theP).id() );
  }

  return theNLightProtonVec;

}

//////////////////////////////////////////////////////////////////////

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

  // NOTE: In the following we account for the full scale restoration 
  // if \mu of the OLP differs from \mu_R.
  // Note: In the GoSam OLP interface, it is possible to directly set 
  // \mu = \mu_R, via the switch SetMuToMuR (for debugging purposes).
  if ( !lastBorn()->hasRunningAlphaS() ) {
    Energy2 muR2 = 
      lastBorn()->renormalizationScale()*
      sqr(lastBorn()->renormalizationScaleFactor());
    if ( muR2 != mu2 ) {
      res -=
	betaZero *
	lastBorn()->orderInAlphaS() * log(muR2/mu2) *
	lastBorn()->me2();
    }
  }

  // include the finite renormalization for DR here; ATTENTION this
  // has to be mentioned in the manual!  see hep-ph/9305239 for
  // details; this guarantees an expansion in alpha_s^\bar{MS} when
  // using dimensional reduction
  if ( isDR() && isDRbar() ) {
    res -= (CA/6.)*lastBorn()->orderInAlphaS()*lastBorn()->me2();
  }

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

//////////////////////////////////////////////////////////////////////

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
describeHerwigDipoleIOperator("Herwig::DipoleIOperator", "Herwig.so");

void DipoleIOperator::Init() {

  static ClassDocumentation<DipoleIOperator> documentation
    ("DipoleIOperator");

  DipoleRepository::registerInsertionIOperator<0,DipoleIOperator>("LightIOperator");

}

