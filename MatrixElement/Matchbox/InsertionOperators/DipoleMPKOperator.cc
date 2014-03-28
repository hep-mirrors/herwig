// -*- C++ -*-
//
// DipoleMPKOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleMPKOperator class.
//

#include "DipoleMPKOperator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;
using Constants::pi;

DipoleMPKOperator::DipoleMPKOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    KQuark(-1.0), KGluon(-1.0) {}

DipoleMPKOperator::~DipoleMPKOperator() {}

IBPtr DipoleMPKOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMPKOperator::fullclone() const {
  return new_ptr(*this);
}

void DipoleMPKOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
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
}

double DipoleMPKOperator::gammaSoft() const {
  double res = (1.+log(1.-x))*PDFx(parton);
  if ( z > x )
    res +=
      (PDFxByz(parton) - z*PDFx(parton)) / (z*(1.-z));
  return res;
}

double DipoleMPKOperator::softLogByz(tcPDPtr p) const {
  double res = ( sqr(log(1.-x))/2. - sqr(pi)/6. ) * PDFx(p);
  if ( z > x ) {
    res += (PDFxByz(p) - z*PDFx(p))*log(1.-z)/(z*(1.-z));
    res -= PDFxByz(p)*log(z)/(z*(1.-z));
  }
  return res;
}

double DipoleMPKOperator::softLog(tcPDPtr p) const {
  double res = sqr(log(1.-x)) * PDFx(p) / 2.;
  if ( z > x ) {
    res += (PDFxByz(p) - z*PDFx(p))*log(1.-z)/(z*(1.-z));
  }
  return res;
}

double DipoleMPKOperator::KBarqq() const {
  assert(abs(parton->id()) < 7);
  double res = 
    2.*softLogByz(parton) +
    (sqr(pi) - 5.)*PDFx(parton);
  if ( z > x ) {
    res += PDFxByz(parton)*( (1.-z) - (1.+z)*log((1.-z)/z) ) / z;
  }
  return (res * CF);
}

double DipoleMPKOperator::KTildeqq() const {
  assert(abs(parton->id()) < 7);
  double res =
    2.*CF*softLog(parton) - CF*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res -= ( CF * (1.+z) * log(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

double DipoleMPKOperator::Pqq() const {
  assert(abs(parton->id()) < 7);
  double res = (3./2.+2.*log(1.-x)) * PDFx(parton);
  if ( z > x ) {
    res += 2.*(PDFxByz(parton) - z*PDFx(parton))/(z*(1.-z));
    res -= PDFxByz(parton) * (1.+z)/z;
  }
  return (CF*res);
}

double DipoleMPKOperator::KBarqg() const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double res = 0.0;
  double factor = CF * ( ( (1.+sqr(1.-z)) / z ) * log((1.-z)/z) + z ) / z;
  int nl= lastBorn()->nLight();
  for ( int f = -lastBorn()->nLight(); f <= nl; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f))*factor;
  }
  return res;
}

double DipoleMPKOperator::KTildeqg() const {
  if ( z < x )
    return 0.0;
  return Pqg() * log(1.-z);
}

double DipoleMPKOperator::Pqg() const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double res = 0.0;
  double factor = CF * ( 1. + sqr(1.-z) ) / sqr(z);
  int nl = lastBorn()->nLight();
  for ( int f = -lastBorn()->nLight(); f <= nl; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f))*factor;
  }
  return res;
}

double DipoleMPKOperator::KBargq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return
    PDFxByz(getParticleData(ParticleID::g)) *
    ( 0.5*(sqr(z)+sqr(1.-z))*log((1.-z)/z) + z*(1.-z) ) / z;
}

double DipoleMPKOperator::KTildegq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return Pgq() * log(1.-z);
}

double DipoleMPKOperator::Pgq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return 0.5 * ( sqr(z) + sqr(1.-z) ) * PDFxByz(getParticleData(ParticleID::g)) / z;
}

double DipoleMPKOperator::KBargg() const {
  assert(parton->id() == ParticleID::g);
  double res = 
    2.* CA* softLogByz(parton) +
    ( CA*( sqr(pi) - 50./9. ) + (8./9.)*lastBorn()->nLight() ) * PDFx(parton);
  if ( z > x ) {
    res += 2.*CA*((1.-z)/z-1.+z*(1.-z))*log((1.-z)/z)*PDFxByz(parton)/z;
  }
  return res;
}

double DipoleMPKOperator::KTildegg() const {
  assert(parton->id() == ParticleID::g);
  double res =
    2.*CA*softLog(parton) - CA*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res += ( 2.*CA * ( (1.-z)/z -1. + z*(1.-z) ) * log(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

double DipoleMPKOperator::Pgg() const {
  assert(parton->id() == ParticleID::g);
  double res = 
    ( (11./6.) * CA - (1./3.) * lastBorn()->nLight() + 2.*CA*log(1.-x) ) * PDFx(parton);
  if ( z > x ) {
    res += 2. * CA * ( PDFxByz(parton) - z*PDFx(parton) ) / (z*(1.-z));
    res += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

//////////////////////////////////////////////////

/**
 * [J^a_{gQ}(z,\mu_Q^2)]_+
 */
double DipoleMPKOperator::Ja_gQplus(double muQ2) const {
  double res = 
    -1. * PDFx(parton) * ( (1.-z)/(2.*(1.-z-muQ2)*(1.-z-muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(1.+muQ2))) );
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * ( (1.-z)/(2.*(1.-z-muQ2)*(1.-z-muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(2.-z+muQ2))) );
  }
  return res;
}

/**
 * [1/(1-z)]_+ * log( (2-z)/(2-z+\mu_Q^2) )
 */
double DipoleMPKOperator::gammaSoft2(double muQ2) const {
  double res = 
    (1./(1.-z)) * ( -1. * PDFx(parton) * 2. * log( 1./(1.+muQ2) ) );
  if ( z > x ) {
    res += (1./(1.-z)) * ( 1./z * PDFxByz(parton) * 2. * log( (2.-z)/(2.-z+muQ2) ) );
  }
  return res;
}

/**
 * The Kscript^{qq}_q contribution
 */
double DipoleMPKOperator::Kscriptqq_q(Energy2 sja, Energy2 mj2) const {
  double muQ2 = mj2/sja;
  double res = 
    2. * softLog(parton) + 
    Ja_gQplus(muQ2) + 
    2. * gammaSoft2(muQ2) - 
    gammaQuark/CF * PDFx(parton) + 
    ( muQ2*log(muQ2/(1.+muQ2)) + 1./2.*(muQ2/(1.+muQ2)) ) * PDFx(parton);
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * (
      -2. * log(2.-z)/(1.-z) );
  }
  return res;
}

/**
 * The Kscript^{qg}_q contribution
 */
double DipoleMPKOperator::Kscriptqg_q(Energy2 sja, Energy2 mj2) const {
  double muQ2 = mj2/sja;
  double res = 0.0;
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * (
      2.*CF/CA*muQ2/z*log(muQ2/(1.-z+muQ2)) );
  }
  return res;
}

/**
 * The Kscript^{gq}_q contribution
 */
double DipoleMPKOperator::Kscriptgq_q(Energy2 sja, Energy2 mj2) const {
  double res = 0.0;
  return res;
}

/**
 * The Kscript^{gg}_q contribution
 */
double DipoleMPKOperator::Kscriptgg_q(Energy2 sja, Energy2 mj2) const {
  double res = Kscriptqq_q(sja,mj2) + CA/CF*Kscriptqg_q(sja,mj2);
  return res;
}

/**
 * J^{a;NS}_{Q\bar{Q}}(\mu_Q^2)
 * Does not include the folding with 1/z*PDF(x/z)*\Theta(z-x)
 */
double DipoleMPKOperator::JaNS_QQ(double muQ2) const {
  double res = 
    10./9. * ( 1. - sqrt(1.-4.*muQ2) ) - 
    8./9. * muQ2 * sqrt(1.-4.*muQ2) + 
    4./3. * log( (1.+sqrt(1.-4.*muQ2))/2. );
  return res;
}

/**
 * [J^a_{Q\bar{Q}}(z,\mu_Q^2)]_{z_+}
 */
double DipoleMPKOperator::Ja_QQzplus(double muQ2) const {
  double zplus = 1. - 4.*muQ2;
  double res = 0.0;
  if ( z > x && z < zplus ) {
    res += ( 1./z*PDFxByz(parton) - 1./zplus*PDFxByzplus(parton,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( z > x && z < zplus ) {
    res += ( 1./zplus*PDFxByzplus(parton,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( zplus > x && z < zplus ) {
    res -= ( 1./zplus*PDFxByzplus(parton,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  return res;
}

/**
 * The Kscript^{qq}_g contribution
 */
double DipoleMPKOperator::Kscriptqq_g(Energy2 sja, Energy2 mj2) const {
  int NLight = lastBorn()->nLight();
  double res = -1.*gammaGluon/CA*gammaSoft();
  for( int f=1; f<=NLight; ++f ) {
    Energy2 mF2 = sqr( getParticleData(f)->mass() );
    double muQ2 = mF2/sja;
    double zplus = 1. - 4*muQ2;
    if( mF2 == ZERO ) continue; // only heavy quarks
    if( sja <= 4.*mF2 ) continue; // sum only over quarks which meet special condition
    res += 1./(2.*CA) * (
      PDFx(parton)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) -
      1./zplus*PDFxByzplus(parton,muQ2)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) + 
      Ja_QQzplus(muQ2) + 
      PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) );
  }
  return res;
}

/**
 * The Kscript^{qg}_g contribution
 */
double DipoleMPKOperator::Kscriptqg_g(Energy2 sja, Energy2 mj2) const {
  double res = 0.0;
  return res;
}

/**
 * The Kscript^{gq}_g contribution
 */
double DipoleMPKOperator::Kscriptgq_g(Energy2 sja, Energy2 mj2) const {
  double res = 0.0;
  return res;
}

/**
 * The Kscript^{gg}_g contribution
 * equals the Kscript^{qq}_g contribution
 */
double DipoleMPKOperator::Kscriptgg_g(Energy2 sja, Energy2 mj2) const {
  int NLight = lastBorn()->nLight();
  double res = -1.*gammaGluon/CA*gammaSoft();
  for( int f=1; f<=NLight; ++f ) {
    Energy2 mF2 = sqr( getParticleData(f)->mass() );
    double muQ2 = mF2/sja;
    double zplus = 1. - 4*muQ2;
    if( mF2 == ZERO ) continue; // only heavy quarks
    if( sja <= 4.*mF2 ) continue; // sum only over quarks which meet special condition
    res += 1./(2.*CA) * (
      PDFx(parton)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) -
      1./zplus*PDFxByzplus(parton,muQ2)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) + 
      Ja_QQzplus(muQ2) + 
      PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) );
  }
  return res;
}

//////////////////////////////////////////////////

double DipoleMPKOperator::PDFx(tcPDPtr pd) const {
  map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator
    cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = make_pair(0.0,0.0);
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  if ( cached->second.first == 0.0 )
    cached->second.first = 
      pdf->xfx(particle,pd,scale,x)/x;
  return cached->second.first;
}

double DipoleMPKOperator::PDFxByz(tcPDPtr pd) const {
  map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator
    cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = make_pair(0.0,0.0);
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  if ( cached->second.second == 0.0 )
    cached->second.second = 
      pdf->xfx(particle,pd,scale,x/z)*z/x;
  return cached->second.second;
}

//////////////////////////////////////////////////

// ! Unsure whether this is actually correct !  //

double DipoleMPKOperator::PDFxByzplus(tcPDPtr pd, double muQ2) const {
  double zplus = 1. - 4.*muQ2;
  map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator
    cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = make_pair(0.0,0.0);
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  if ( cached->second.second == 0.0 )
    cached->second.second = 
      pdf->xfx(particle,pd,scale,x/zplus)*zplus/x;
  return cached->second.second;
}

//////////////////////////////////////////////////

bool DipoleMPKOperator::apply(tcPDPtr pd) const {
  // DipoleMPKOperator has to apply also to massive partons
  return
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

bool DipoleMPKOperator::apply(const cPDVector& pd) const {
  if ( !apply(pd[0]) && !apply(pd[1]) )
    return false;
  bool first = false;
  bool second = false;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {

    if ( !first ) {
      if ( apply(*p) )
	first = true;
    }
    else {
      if ( apply(*p) )
	second = true;
    }

  }
  return first && second;
}

double DipoleMPKOperator::sumParton(int id) const {

// sumParton(int id) sums for a specific id=a'=0,1
// over all a

  pdf =
    id == 0 ?
    lastXCombPtr()->partonBins().first->pdf() :
    lastXCombPtr()->partonBins().second->pdf();

  x = 
    id == 0 ?
    lastXCombPtr()->lastX1() :
    lastXCombPtr()->lastX2();

  parton =
    id == 0 ?
    mePartonData()[0] :
    mePartonData()[1];

  particle =
    id == 0 ?
    lastParticles().first->dataPtr() :
    lastParticles().second->dataPtr();

  using namespace RandomHelpers;

  double r = insertionRandomNumbers().front();
  double eps = 1e-3;

  pair<double,double> zw =
    generate((piecewise(),
	      flat(0.0,x),
	      match(inverse(0.0,x,1.0) + inverse(1.0+eps,x,1.0))),r);

  z = zw.first;
  double mapz = zw.second;

  for ( map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator cache =
	  pdfCache.begin(); cache != pdfCache.end(); ++cache )
    cache->second = make_pair(0.0,0.0);

  assert(pdf);

  double res = 0.;

  // The apply function makes really only a difference in the sum over j
  // (sum over finals):
  // It allows for massive particles in general, but we trust that the
  // incoming partons are massless. Those contributions that aren't inside
  // the sum over j are the same in the massless and the massive P+K, and
  // should not be accounted here at all.
  
  // We could have a g(a)->Q(a')Qbar(j) or a g(a)->Qbar(a')Q(j) splitting
  // in the initial state. This should then also be correctly handled.

  ////////////////////////////
  // K operator             //
  // non-color correlated   //
  ////////////////////////////

//   if ( mePartonData()[id]->id() == ParticleID::g )
//     res += (KBargg() + KBarqg())*lastBorn()->me2();
// 
//   if ( abs(mePartonData()[id]->id()) < 7 )
//     res += (KBarqq() + KBargq())*lastBorn()->me2();
// 
//   //////////
// 
//     // The non-color correlated contributions of K
//     // are equal in the massive and massless case.
// 
//     double zeroMassKNonCorrelated = res;
// 
//   //////////

  ////////////////////////////

  double theGammaSoft = 0.0;

  double thePqq = 0.0;
  double thePqg = 0.0;
  double thePgq = 0.0;
  double thePgg = 0.0;

  double ifCorrelated = 0.0;
  double fiCorrelated = 0.0;

  int idi = 2;
  vector<Lorentz5Momentum>::const_iterator Pi = meMomenta().begin() + 2;

  double disFinite = 0.0;
  double glueFinite = 0.0;

  //////////////////////////////////////////////////////////
  // Initial-Final and Final-Initial contributions        //
  //////////////////////////////////////////////////////////

  //////////

    // Prepare also to cache the total contribution to K and
    // P of the mi=0 case, after the sum over i.

    double zeroMassK = 0.0;
    double zeroMassP = 0.0;

  //////////

  for ( cPDVector::const_iterator i = mePartonData().begin() + 2;
	i != mePartonData().end(); ++i, ++Pi, ++idi ) {

    if ( !apply(*i) || lastBorn()->noDipole(idi,id) )
      continue;

    fiCorrelated = lastBorn()->colourCorrelatedME2(make_pair(idi,id));
    ifCorrelated = lastBorn()->colourCorrelatedME2(make_pair(id,idi));

    if ( theGammaSoft == 0.0 )
      theGammaSoft = gammaSoft();

    if ( mePartonData()[id]->id() == ParticleID::g &&
	 thePgg == 0.0 ) {
      thePqg = Pqg();
      thePgg = Pgg();
    }

    if ( abs(mePartonData()[id]->id()) < 7 &&
	 thePqq == 0.0 ) {
      thePgq = Pgq();
      thePqq = Pqq();
    }

    ////////////////
    // K operator //
    ////////////////

    //////////

      // Last term in massless K operator in (C.31) in massless paper
      zeroMassK +=
        ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
        theGammaSoft * fiCorrelated;

    //////////

    Energy2 mj2 = sqr((**i).mass());
    Energy2 sja = 2.*( (*Pi) * meMomenta()[id]/z );

    // Next-to-Last term in massive K operator in (6.55) in massive paper
    // corresponds in the massless limit to the last term of the massless 
    // K operator in (C.31) in massless paper
    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res -=
        ( (**i).id() == ParticleID::g ? CA : CF ) *
        fiCorrelated * 
        ( (**i).id() == ParticleID::g ? ( Kscriptqg_g(sja,mj2) + Kscriptgg_g(sja,mj2) ) : ( Kscriptqg_q(sja,mj2) + Kscriptgg_q(sja,mj2) ) );
    }
    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res -=
        ( (**i).id() == ParticleID::g ? CA : CF ) *
        fiCorrelated * 
        ( (**i).id() == ParticleID::g ? ( Kscriptqq_g(sja,mj2) + Kscriptgq_g(sja,mj2) ) : ( Kscriptqq_q(sja,mj2) + Kscriptgq_q(sja,mj2) ) );
    }

    // The regular splitting functions 
    // without the folding with 1/z*PDF(x/z)*\Theta(z-x)
    double thePqqreg = -1.*CF*(1.+z);
    double thePqgreg = CF*(1.+(1.-z)*(1.-z))/z;
    double thePgqreg = 1./2.*(z*z + (1.-z)*(1.-z));
    double thePggreg = 2.*CA*((1.-z)/z - 1. + z*(1.-z));

    // Last term of the massive K operator in (6.55) massive paper
    // Vanishes in the massless limit
    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res -=
        ifCorrelated * 
        ( ( z>x ? 1./z*PDFxByz(parton)*(thePqgreg+thePggreg)*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
          (gammaGluon)*PDFx(parton)*( log((sja-2.*sqrt(mj2)*sqrt(sja+mj2)+2.*mj2)/sja) + 2.*sqrt(mj2)/(sqrt(sja+mj2)+sqrt(mj2)) ) );
    }
    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res -=
        ifCorrelated * 
        ( ( z>x ? 1./z*PDFxByz(parton)*(thePqqreg+thePgqreg)*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
          (gammaQuark)*PDFx(parton)*( log((sja-2.*sqrt(mj2)*sqrt(sja+mj2)+2.*mj2)/sja) + 2.*sqrt(mj2)/(sqrt(sja+mj2)+sqrt(mj2)) ) );
    }

    ////////////////
    // P operator //
    ////////////////

    // P(m) = P(m=0).

    // Still need to regulate that the DipoleIOperator
    // contains the integrated finite extra terms that
    // render the massless Dipoles positive definite.
    // At the moment the same extra terms are also in
    // use in the massive dipoles. This has yet to be 
    // clarified, though, for the massive case.

    double theLog = log(scale/(2.*((*Pi)*meMomenta()[id])));

    //////////

      if ( mePartonData()[id]->id() == ParticleID::g ) {
        zeroMassP +=
          ( thePgg + thePqg ) * theLog * ifCorrelated;
      }

      if ( abs(mePartonData()[id]->id()) < 7 ) {
        zeroMassP +=
          ( thePqq + thePgq ) * theLog * ifCorrelated;
        if ( disFinite == 0.0 && z > x ) {
          disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
        }
        if ( z > x )
          zeroMassP -= disFinite*ifCorrelated;
      }

      if ( abs(mePartonData()[idi]->id()) < 7 ) {
        if ( disFinite == 0.0 && z > x ) {
          disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
        }
        if ( z > x )
          zeroMassP -= disFinite*fiCorrelated;
      }

      if ( mePartonData()[idi]->id() == ParticleID::g ) {
        if ( glueFinite == 0.0 && z > x ) {
          glueFinite = 2.*CA*PDFxByz(parton)*(1.+z/6.)/z;
        }
        if ( z > x )
          zeroMassP -= glueFinite*fiCorrelated;
      }

    //////////

    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res +=
	( thePgg + thePqg ) * theLog * ifCorrelated;
    }

    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res +=
	( thePqq + thePgq ) * theLog * ifCorrelated;
      if ( disFinite == 0.0 && z > x ) {
	disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      if ( z > x )
	res -= disFinite*ifCorrelated;
    }

    if ( abs(mePartonData()[idi]->id()) < 7 ) {
      if ( disFinite == 0.0 && z > x ) {
	disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      if ( z > x )
	res -= disFinite*fiCorrelated;
    }

    if ( mePartonData()[idi]->id() == ParticleID::g ) {
      if ( glueFinite == 0.0 && z > x ) {
	glueFinite = 2.*CA*PDFxByz(parton)*(1.+z/6.)/z;
      }
      if ( z > x )
	res -= glueFinite*fiCorrelated;
    }

  } // end loop over i

  //////////////////////////////////////////////////////////
  // Initial-Initial contributions                        //
  //////////////////////////////////////////////////////////

//   // The Initial-Initial contributions of P+K are the same
//   // in the massive and massless case.
// 
//   //////////
// 
//     double zeroMassII = 0.0;
// 
//   //////////
// 
//   if ( mePartonData()[ id == 0 ? 1 : 0 ]->coloured() &&
//        !lastBorn()->noDipole(id == 0 ? 0 : 1,
// 			     id == 0 ? 1 : 0) ) {
// 
//     if ( mePartonData()[id]->id() == ParticleID::g &&
// 	 thePgg == 0.0 ) {
//       thePqg = Pqg();
//       thePgg = Pgg();
//     }
// 
//     if ( abs(mePartonData()[id]->id()) < 7 &&
// 	 thePqq == 0.0 ) {
//       thePgq = Pgq();
//       thePqq = Pqq();
//     }
// 
//     double theLog = log(scale/(2.*(meMomenta()[0]*meMomenta()[1])));
// 
//     pair<int,int> corr = id == 0 ? make_pair(0,1) : make_pair(1,0);
//     double iiCorrelated = lastBorn()->colourCorrelatedME2(corr);
// 
//     //////////
// 
//       if ( mePartonData()[id]->id() == ParticleID::g ) {
//         zeroMassII +=
//           ( thePgg + thePqg ) * theLog * iiCorrelated;
//         zeroMassII -=
//           ( KTildegg() + KTildeqg() ) * iiCorrelated;
//       }    
// 
//       if ( abs(mePartonData()[id]->id()) < 7 ) {
//         zeroMassII +=
//           ( thePqq + thePgq ) * theLog * iiCorrelated;
//         zeroMassII -=
//           ( KTildeqq() + KTildegq() ) * iiCorrelated;
//       }
// 
//     //////////
// 
//     if ( mePartonData()[id]->id() == ParticleID::g ) {
//       res +=
// 	( thePgg + thePqg ) * theLog * iiCorrelated;
//       res -=
// 	( KTildegg() + KTildeqg() ) * iiCorrelated;
//     }    
// 
//     if ( abs(mePartonData()[id]->id()) < 7 ) {
//       res +=
// 	( thePqq + thePgq ) * theLog * iiCorrelated;
//       res -=
// 	( KTildeqq() + KTildegq() ) * iiCorrelated;
//     }
// 
//   }

  //////////////////////////////////////////////////////////

  // Subtract massless contributions if needed

  idi = 2;
  for ( cPDVector::const_iterator i = mePartonData().begin() + 2;
	i != mePartonData().end(); ++i, ++idi ) {
    if (mePartonData()[idi]->mass() == ZERO) {
//       res -= ( zeroMassII + 
//                zeroMassKNonCorrelated + 
//                zeroMassK + 
//                zeroMassP );
      res -= ( zeroMassK + 
               zeroMassP );
    }
  }

  //////////////////////////////////////////////////////////

  return res * mapz;

}

double DipoleMPKOperator::me2() const {

  scale = lastBorn()->lastScale();

  double res = 0.0;

  if ( apply(mePartonData()[0]) ) {
    if ( mePartonData()[0]->coloured() ) {
      if ( mePartonData()[1]->coloured() )
	res += lastBorn()->pdf2()*sumParton(0);
      else
	res += sumParton(0);
    }
  }

  if ( apply(mePartonData()[1]) ) {
    if ( mePartonData()[1]->coloured() ) {
      if ( mePartonData()[0]->coloured() )
	res += lastBorn()->pdf1()*sumParton(1);
      else
	res += sumParton(1);
    }
  }

  return (lastBorn()->lastAlphaS()/(2.*pi)) * res;

}

void DipoleMPKOperator::persistentOutput(PersistentOStream & os) const {
  os << CA << CF << gammaQuark << gammaGluon << KQuark << KGluon
     << ounit(scale,GeV2) << pdf << particle << x << z << pdfCache << parton;
}

void DipoleMPKOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> KQuark >> KGluon
     >> iunit(scale,GeV2) >> pdf >> particle >> x >> z >> pdfCache >> parton;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipoleMPKOperator,MatchboxInsertionOperator>
describeHerwigDipoleMPKOperator("Herwig::DipoleMPKOperator", "Herwig.so");

void DipoleMPKOperator::Init() {

  static ClassDocumentation<DipoleMPKOperator> documentation
    ("DipoleMPKOperator");

  DipoleRepository::registerInsertionOperator<0,DipoleMPKOperator>("MassivePKOperator");

}

