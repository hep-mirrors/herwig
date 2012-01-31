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
    KQuark(-1.0), KGluon(-1.0),
    theUseDR(false) {}

DipoleMIOperator::~DipoleMIOperator() {}

IBPtr DipoleMIOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMIOperator::fullclone() const {
  return new_ptr(*this);
}

void DipoleMIOperator::dumpInfo(const string& prefix) const {
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
  generator()->log() << "  | Born ME\n";
  lastBorn()->dumpInfo(prefix+"  | ");
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

  bool oneMassive = false;
  bool first = false;
  bool second = false;
  int i=0;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p, ++i ) {
    if ( abs((**p).id()) < 7&& (**p).mass() != ZERO )
      oneMassive = true;
    if ( !first ) {
      if ( apply(*p, i<2 ? true : false) )
	first = true;
    } else {
      if ( apply(*p, i<2 ? true : false) )
	second = true;
    }
  }
  return first && second && oneMassive;
}

bool DipoleMIOperator::apply(tcPDPtr pd, bool incoming) const {
  return
//     pd->mass() == ZERO &&
//     (abs(pd->id()) < 6 || pd->id() == ParticleID::g);
    ( abs(pd->id()) < 6 && !incoming && pd->mass() != ZERO ) ||
    ( abs(pd->id()) < 6 && incoming && pd->mass() == ZERO ) ||
    pd->id() == ParticleID::g;
}

double DipoleMIOperator::me2() const {

  double res = 0.;

  int idj = 0; int idk = 0;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j, ++idj ) {

    if ( !apply(*j) )
      continue;

    idk = 0;

    for ( cPDVector::const_iterator k = mePartonData().begin();
	  k != mePartonData().end(); ++k, ++idk ) {

      if ( !apply(*k) )
	continue;

      if ( j == k || lastBorn()->noDipole(idj,idk) )
	continue;
      
      // new
      Energy2 sjk = 2.*meMomenta()[idj]*meMomenta()[idk];
      double kappa=0.;
      // renormalization scale (should cancel with something)
      Energy2 mu2 = lastBorn()->renormalizationScale();
      
      // ASSUMING INCOMING PARTICLES HAVE ARE AT 0 AND 1 IN mePartonData, OUTGOING >= 2
      // this is the I operator in the case of no initial-state hadrons (6.16)
      if( idj>=2 && idk>=2 )
	res +=
	  ( ((**j).id() == ParticleID::g ? CA : CF) * (Vj(**j,**k,sjk,kappa)-sqr(pi)/3.) +
	    ((**j).id() == ParticleID::g ? GammaGluon() : GammaQuark(**j)) +
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * (1.+log(mu2/sjk)) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) ) *
	  lastBorn()->colourCorrelatedME2(make_pair(idj,idk));
      // additional terms in the case of one initial-state hadron (6.51)
      if( idj>=2 && idk<2 ) {
	res +=
	  ( ((**j).id() == ParticleID::g ? CA : CF) * (Vj(**j,**k,sjk,kappa)-sqr(pi)/3.) +
	    ((**j).id() == ParticleID::g ? GammaGluon() : GammaQuark(**j)) +
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * (1.+log(mu2/sjk)) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) ) *
	  lastBorn()->colourCorrelatedME2(make_pair(idj,idk));
	res +=
	  ( ((**k).id() == ParticleID::g ? CA : CF) * (Vj(**k,**j,sjk,2./3.,true)-sqr(pi)/3.) +
	    // MAKE SURE lastScale() RETURNS (pi+pj)^2, NOT  2.*pi*pj !!!!
	    ((**k).id() == ParticleID::g ? gammaGluon : gammaQuark) * ( log(lastScale()/mu2) + log(mu2/sjk) + 1. ) +
	    ((**k).id() == ParticleID::g ? KGluon : KQuark) ) *
	  lastBorn()->colourCorrelatedME2(make_pair(idk,idj));
      }
      // additional terms in the case of two initial-state hadrons (6.66)
      if( idj<2 && idk<2 )
	res +=
	  ( ((**j).id() == ParticleID::g ? CA : CF) * ( 0.5*sqr(log(lastScale()/sjk)) - sqr(pi)/3. ) +
	    ((**j).id() == ParticleID::g ? gammaGluon : gammaQuark) * ( log(lastScale()/sjk) + 1. ) +
	    ((**j).id() == ParticleID::g ? KGluon : KQuark) ) *
	  lastBorn()->colourCorrelatedME2(make_pair(idj,idk));
      
    }
  }

  res *= ( - lastBorn()->lastAlphaS() / (2.*pi) );
  
  if( gsl_isnan(res) ) cout << "DipoleMIOperator::me2 nan" << endl;
  if( gsl_isinf(res) ) cout << "DipoleMIOperator::me2 inf" << endl;
  
  // MAKE SURE lastScale() RETURNS (pi+pj)^2, NOT  2.*pi*pj !!!!
  // check against e+ e- --> q qbar
  
//   Energy2 mj2 =  sqr(mePartonData()[appliedParton]->mass());
//   Energy2 sjk = 2.*meMomenta()[2]*meMomenta()[3];
//   Energy2 Qjk2 = sjk + 2.*mj2;
//   
//   double v = sqrt(1.-4.*mj2/Qjk2);
//   double VS = (1.+v*v)/(2.*v)*log(2./(1.+v*v))*log((1.-v)/(1.+v)) - 
//     (1.+v*v)/(4.*v)*sqr(log((1.-v)/(1.+v))) - (1.+v*v)/(2.*v)*sqr(pi)/6. +
//     (1.+v*v)/(2.*v)*log(2./(1.+v*v))*log((1.-v)/(1.+v));
//   double VNS = 3./2.*log((1.+v*v)/2.) +
//     (1.+v*v)/(2.*v) * ( 2.*log((1.-v)/(1.+v))*log(2.*(1.+v*v)/sqr(1.+v)) +
//       2.*gsl_sf_dilog(sqr((1.-v)/(1.+v))) - 2.*gsl_sf_dilog(2.*v/(1.+v)) - sqr(pi)/6. ) +
//     log(1.-1./2.*sqrt(1.-v*v)) - 2.*log(1.-sqrt(1.-v*v)) -
//       (1.-v*v)/(1.+v*v)*log(sqrt(1.-v*v)/(2.-sqrt(1.-v*v))) -
//     sqrt(1.-v*v)/(2.-sqrt(1.-v*v)) + 2.*(1.-v*v-sqrt(1.-v*v))/(1.+v*v) + sqr(pi)/2.;
//     
//     
//   double target = 2.*lastBorn()->lastAlphaS() / (2.*pi) * CF *
//     ( VS + VNS + 1./2.*log((1.-v*v)/4.) - 2. + 3./2.*log(2./(1+v*v)) -
//       sqr(pi)/3. + (gammaQuark+KQuark)/CF );
//   
//   target *= lastBorn()->me2();
//    
//   cout << "+++++++++++++++++++++++++++++MI diff resc " << (res-target)/(res+target) << endl;
//   cout << "res " << res << "  target " << target << " @id " << mePartonData()[appliedParton]->id() << endl;
//   cout << "VS+VNS part " << 2.*lastBorn()->lastAlphaS() / (2.*pi) * CF * ( VS + VNS ) <<
//     "  while VS+VNS " << VS+VNS << endl;
  
  return res;

}

double DipoleMIOperator::GammaQuark(const ParticleData& j) const {
  // renormalization scale (should cancel with something)
  Energy2 mu2 = lastBorn()->renormalizationScale();
  // massless quark: only finite remainder
  // MAKE SURE lastScale() RETURNS (pi+pj)^2, NOT  2.*pi*pj !!!!
  if( j.mass() == ZERO ) return gammaQuark * log(lastScale()/mu2);
  // massive quark
  return CF * ( log(lastScale()/mu2) + 0.5*log(sqr(j.mass())/mu2) - 2. );
}

double DipoleMIOperator::GammaGluon() const {
  // main contribution cancels with VjNS, only finite 1/eps remainder
  return gammaGluon * log( lastScale() / lastBorn()->renormalizationScale() );
}

double DipoleMIOperator::Vj(const ParticleData& j, const ParticleData& k,
			    Energy2 sjk, double kappa, bool mFSetEmpty) const {
  
  double res = 0.;
  
  Energy2 mj2 = sqr(j.mass()), mk2 = sqr(k.mass());
//   Energy2 sjk = 2.*meMomenta()[idj]*meMomenta()[idk];
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
  if( mj2 == ZERO && mk2 == ZERO )
    // only finite 1/eps^2 remainder
    res += 1./2.*sqr(log(Qjk2/sjk));
  // one mass zero
  else if( mj2 == ZERO || mk2 == ZERO ) {
    Energy2 m2 = sqr(l.mass());
    res += -1./4.*sqr(log(m2/sjk)) - sqr(pi)/12. -
      1./2.*log(m2/sjk)*log(sjk/Qjk2) - 1./2.*log(m2/Qjk2)*log(sjk/Qjk2);
    // finite 1/eps^2 remainder
    res += 1./4.*sqr(log(Qjk2/sjk));
    // finite 1/eps remainder
    res += 1./2.*log(m2/sjk)*log(Qjk2/sjk);
  }
  // no mass zero
  else if( mj2 != ZERO && mk2 != ZERO ) {
    res += 1./vjk * ( -1./4.*sqr(log(rhoj*rhoj)) - 1./4.*sqr(log(rhok*rhok)) -
      sqr(pi)/6. + ( rho==0. ? 0. : log(rho)*log(Qjk2/sjk) ) );
    // finite 1/eps remainder
    res += 1./vjk * ( rho==0. ? 0. : log(rho)*log(Qjk2/sjk) );
  }
  else cout << "problem occurred in DipoleMIOperator::Vj -- S part" << endl;
  
//   double resS = res;
  
  if(gsl_isnan(res)) cout << "Vj S nan" << "  j " << j.id() << "  k " << k.id() << endl;
  
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
	  for( int f=1; f< 7; ++f ) {
	    Energy2 mF2 = sqr( getParticleData(f)->mass() );
	    // only heavy quarks
	    if( mF2 == ZERO ) continue;
	    // sum only over quarks which meet special condition
	    if( sjk <= 4.*sqrt(mF2)*(sqrt(mF2)+mk) )
	      continue;
	    double rho1 = sqrt( 1. - 4.*mF2 / sqr(Qjk-mk) );
// 	    double rho2 = sqrt( 1. - 4.*mF2 / (Qjk2-mk2) );
	    res += 2./3./CA * ( log((1.+rho1)/2.) - rho1/3.*(3.+sqr(rho1)) - 0.5*log(mF2/Qjk2) );
	  }
      }
    }
    // k is massive quark
    else {
      assert( abs(k.id()) < 7);
      // part common to j massless quark or gluon
//     if( mk > ZERO )
      res += sqr(pi)/6. - gsl_sf_dilog(sjk/Qjk2);
      // j is massless (incoming) quark
      if( abs(j.id()) < 7)
	res += gammaQuark/CF * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) );
      // j is gluon
      else if( j.id() == ParticleID::g ) {
	// part independent of other heavy quark flavours
  //       if( mk > ZERO )
	res += gammaGluon/CA * ( log(sjk/Qjk2) - 2.*log((Qjk-mk)/Qjk) - 2.*mk/(Qjk+mk) ) +
	  (kappa-2./3.) * mk2/sjk * (1./CA*lastBorn()->nLight()-1.) * log(2.*mk/(Qjk+mk));
	// part containing other heavy quark flavours
	if( !mFSetEmpty )
	  for( int f=1; f< 7; ++f ) {
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
  os << CA << CF << gammaQuark << gammaGluon << KQuark << KGluon << theUseDR;
}

void DipoleMIOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> KQuark >> KGluon >> theUseDR;
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

