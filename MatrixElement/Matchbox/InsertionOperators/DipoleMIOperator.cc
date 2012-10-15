// -*- C++ -*-
//
// DipoleMIOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMIOperator class.
//

#include "DipoleMIOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"

#include <gsl/gsl_sf_dilog.h>
// TODO: remove
// only for checking for NaN or inf
#include <gsl/gsl_math.h>

using namespace Herwig;
using Constants::pi;

DipoleMIOperator::DipoleMIOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    KQuark(-1.0), KGluon(-1.0) {}

DipoleMIOperator::~DipoleMIOperator() {}

IBPtr DipoleMIOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMIOperator::fullclone() const {
  return new_ptr(*this);
}

void DipoleMIOperator::setBorn(Ptr<MatchboxMEBase>::tptr me) {
  MatchboxInsertionOperator::setBorn(me);
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
  gammaQuark = (3./2.)*CF;
  gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLight();
  KQuark = (7./2.-sqr(pi)/6.)*CF;
  KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLight();
  if ( isDR() ) {
    gammaQuark -= CF/2.;
    gammaGluon -= CA/6.;
  }
}

bool DipoleMIOperator::apply(const cPDVector& pd) const {
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {
    for ( cPDVector::const_iterator q = pd.begin();
	  q != pd.end(); ++q ) {
      if ( p == q )
	continue;
      if ( apply(*p,*q) )
	return true;
    }
  }
  return false;
}

bool DipoleMIOperator::apply(tcPDPtr pd1, tcPDPtr pd2) const {
  return
    (abs(pd1->id()) < 7 || pd1->id() == ParticleID::g) &&
    (abs(pd2->id()) < 7 || pd2->id() == ParticleID::g) &&
    (pd1->mass() != ZERO || pd2->mass() != ZERO);
}

double DipoleMIOperator::me2() const {

  if ( !isCS() )
    throw InitException() <<
      "DipoleMIOperator only implemented in the Catani-Seymour scheme.";

  double res = 0.;

  int idj = 0; int idk = 0;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j, ++idj ) {

    idk = 0;

    for ( cPDVector::const_iterator k = mePartonData().begin();
	  k != mePartonData().end(); ++k, ++idk ) {

      if ( j == k || lastBorn()->noDipole(idj,idk) )
	continue;

      if ( !apply(*j,*k) )
	continue;

      // NOTE: massless dipoles handled by DipoleIOperator
      // sum over massive flavours occurs in j==gluon contribution
      if ( abs((*j)->id()) < 6 && (*j)->mass() == ZERO &&
	   (*k)->mass() == ZERO )
 	continue;
      // NOTE: j,k = incoming same contribution as in DipoleIOperator
      if ( idj < 2 && idk < 2 )
	continue;
      
      Energy2 sjk = 2.*meMomenta()[idj]*meMomenta()[idk];
      double kappa=0.;
      
      res +=
	( ((**j).id() == ParticleID::g ? CA : CF) *
	  ( ( idj >= 2 ? Vj(**j,**k,sjk,kappa) : Vj(**j,**k,sjk,2./3.,true) )
	    - sqr(pi)/3. ) +
	  ((**j).id() == ParticleID::g ? GammaGluon() : GammaQuark(**j,sjk)) +
	  // factor (1.+log(mu2/sjk)) absorbed in Gamma_j
	  ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) +
	  ((**j).id() == ParticleID::g ? KGluon : KQuark) ) *
	lastBorn()->colourCorrelatedME2(make_pair(idj,idk));
      // contributions counted here AND in DipoleIOperator if j,k, massless
      // (here only if j==gluon)
      // subtract DipoleIOperator contribution.
      if ( (**j).mass() == ZERO && (**k).mass() == ZERO )
	res -=
	  ( ((**j).id() == ParticleID::g ? CA : CF) * (-sqr(pi)/3.) +
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) ) *
	  lastBorn()->colourCorrelatedME2(make_pair(idj,idk));
      
    }
  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );
  
  if( gsl_isnan(res) ) cout << "DipoleMIOperator::me2 nan" << endl;
  if( gsl_isinf(res) ) cout << "DipoleMIOperator::me2 inf" << endl;
  
  return res;

}

double DipoleMIOperator::GammaQuark(const ParticleData& j, Energy2 sjk) const {
  if ( j.mass() == ZERO )
    return 0.;
  // massive quark, last term see above!
  // CF * (-log(mu2/sjk) + 1/2 log(mj2/mu2) - 2) + 3/2 CF * log(mu2/sjk)
  return CF * ( 0.5*log(sqr(j.mass())/sjk) - 2. );
}

// TODO: kill
double DipoleMIOperator::GammaGluon() const {
  // main contribution cancels with VjNS, only finite 1/eps remainder
  return 0.;
}

// NOTE: no finite remainder of epsilon poles here.
double DipoleMIOperator::Vj(const ParticleData& j, const ParticleData& k,
			    Energy2 sjk, double kappa, bool mFSetEmpty) const {
  
  double res = 0.;
  
  Energy2 mj2 = sqr(j.mass()), mk2 = sqr(k.mass());
  Energy2 Qjk2 = sjk + mj2 + mk2;
  Energy Qjk = sqrt(Qjk2);
  Energy mj = j.mass(), mk = k.mass();
  
  double vjk = rootOfKallen(Qjk2,mj2,mk2) / sjk;
  double rho = sqrt( abs(1.-vjk)/(1.+vjk) ); // abs() because for small mass 1.-vjk can get O(-1.e-16)
  double rhoj = sqrt( ( 1. - vjk + 2.*mj2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) /
    ( 1. + vjk + 2.*mj2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) );
  double rhok = sqrt( ( 1. - vjk + 2.*mk2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) /
    ( 1. + vjk + 2.*mk2/Qjk2 / (1.-mj2/Qjk2-mk2/Qjk2) ) );
  
  ParticleData l = ( mj2 == ZERO ? k : j );
  
  // S part (6.20)
  
  // both masses zero
  if( mj2 == ZERO && mk2 == ZERO );
  // one mass zero
  else if( mj2 == ZERO || mk2 == ZERO ) {
    Energy2 m2 = sqr(l.mass());
    res += -1./4.*sqr(log(m2/sjk)) - sqr(pi)/12. -
      1./2.*log(m2/sjk)*log(sjk/Qjk2) - 1./2.*log(m2/Qjk2)*log(sjk/Qjk2);
  }
  // no mass zero
  else if( mj2 != ZERO && mk2 != ZERO ) {
    res += 1./vjk * ( -1./4.*sqr(log(rhoj*rhoj)) - 1./4.*sqr(log(rhok*rhok)) -
      sqr(pi)/6. + ( rho==0. ? 0. : log(rho)*log(Qjk2/sjk) ) );
  }

  // NS part (6.21)-(6.26)
  
  // V_q (j is quark)
  // j is massive quark
  if( mj2 != ZERO ) {
    assert( abs(j.id()) < 7);
    res += gammaQuark/CF * log(sjk/Qjk2); // iff j and/or k is massive quark (6.21),(6.22)
    // k is massive quark (6.21)
    if( mk2 != ZERO ) {
      assert( abs(k.id()) < 7);
      res += 1./vjk * ( ( rho==0. ? 0. : log(rho*rho)*log(1.+rho*rho) ) + 2.*gsl_sf_dilog(rho*rho) -
	gsl_sf_dilog(1.-rhoj*rhoj) - gsl_sf_dilog(1.-rhok*rhok) - sqr(pi)/6. ) +
	log((Qjk-mk)/Qjk) - 2.*log((sqr(Qjk-mk)-mj2)/Qjk2) - 2.*mj2/sjk*log(mj/(Qjk-mk)) -
	mk/(Qjk-mk) + 2.*mk*(2.*mk-Qjk)/sjk + sqr(pi)/2.;
    }
    // k is massless parton (6.22)
    else {
      res += sqr(pi)/6. - gsl_sf_dilog(sjk/Qjk2) -
	2.*log(sjk/Qjk2) - mj2/sjk*log(mj2/Qjk2);
    }
  }
  // V_q / V_g (j is massless parton)
  else {
    // k is massless parton
    if( mk == ZERO ) {
      // only contribution if j is gluon
      if( j.id() == ParticleID::g ) {
	// sum over all quark flavours
	if( !mFSetEmpty )
	  // TODO: make fmax depend on matrix element
	  for( int f=1; f< 6; ++f ) {
	    Energy2 mF2 = sqr( getParticleData(f)->mass() );
	    // only heavy quarks
	    if( mF2 == ZERO ) continue;
	    // sum only over quarks which meet special condition
	    if( sjk <= 4.*sqrt(mF2)*(sqrt(mF2)+mk) )
	      continue;
	    double rho1 = sqrt( 1. - 4.*mF2 / sqr(Qjk-mk) );
	    res += 2./3./CA * ( log((1.+rho1)/2.) - rho1/3.*(3.+sqr(rho1)) - 0.5*log(mF2/Qjk2) );
	  }
      }
    }
    // k is massive quark
    else {
      assert( abs(k.id()) < 7);
      // part common to j massless quark or gluon
      res += sqr(pi)/6. - gsl_sf_dilog(sjk/Qjk2);
      // j is massless (incoming) quark
      if( abs(j.id()) < 7)
	res += gammaQuark/CF * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) );
      // j is gluon
      else if( j.id() == ParticleID::g ) {
	// part independent of other heavy quark flavours
	res += gammaGluon/CA * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) ) +
	  (kappa-2./3.) * mk2/sjk * (1./CA*lastBorn()->nLight()-1.) * log(2.*mk/(Qjk+mk));
	// part containing other heavy quark flavours
	if( !mFSetEmpty )
	  // TODO: make fmax dependent on matrix element
	  for( int f=1; f< 6; ++f ) {
	    Energy2 mF2 = sqr( getParticleData(f)->mass() );
	    // only heavy quarks
	    if( mF2 == ZERO ) continue;
	    // sum only over quarks which meet special condition
	    if( sjk <= 4.*sqrt(mF2)*(sqrt(mF2)+mk) )
	      continue;
	    double rho1 = sqrt( 1. - 4.*mF2 / sqr(Qjk-mk) );
	    double rho2 = sqrt( 1. - 4.*mF2 / (Qjk2-mk2) );
	    res += 2./3./CA * ( log((Qjk-mk)/Qjk) + mk*rho1*rho1*rho1/(Qjk+mk) + log((1.+rho1)/2.) -
	      rho1/3.*(3.+sqr(rho1)) - 1./2.*log(mF2/Qjk2) ) +
	      1./CA * ( rho2*rho2*rho2*log((rho2-rho1)/(rho2+rho1)) - log((1.-rho1)/(1.+rho1)) -
	      8.*rho1*mF2/sjk );
	  }
      }
    }
  }
    
  return res;
}
  
void DipoleMIOperator::persistentOutput(PersistentOStream & os) const {
  os << CA << CF << gammaQuark << gammaGluon << KQuark << KGluon;
}

void DipoleMIOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> KQuark >> KGluon;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMIOperator,MatchboxInsertionOperator>
describeHerwigDipoleMIOperator("Herwig::DipoleMIOperator", "HwMatchbox.so");

void DipoleMIOperator::Init() {

  static ClassDocumentation<DipoleMIOperator> documentation
    ("DipoleMIOperator");

  DipoleRepository::registerInsertionOperator<DipoleMIOperator>("MassiveIOperator");

}

