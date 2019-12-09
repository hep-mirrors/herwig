// -*- C++ -*-
//
// DipoleMPKOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
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
#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

#include <gsl/gsl_sf_dilog.h>

using namespace Herwig;
using Constants::pi;

DipoleMPKOperator::DipoleMPKOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    KQuark(-1.0), KGluon(-1.0),
    scale(ZERO), x(0.), z(0.) {}

DipoleMPKOperator::~DipoleMPKOperator() {}

IBPtr DipoleMPKOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleMPKOperator::fullclone() const {
  return new_ptr(*this);
}

//////////////////////////////////////////////////////////////////////

bool DipoleMPKOperator::apply(const cPDVector& pd) const {

  // DipoleMPKOperator should apply as soon as massive 
  // partons can occur in the overall process.
  // DipolePKOperator should not apply then.

  if ( !apply(pd[0]) && !apply(pd[1]) ) {
    return false;
  }

  // A gluon in the Born final state can give rise to 
  // a splitting g->Q\bar{Q} in the real radiation.
  // This can happen if massive partons are specified
  // inside the jet (aka if the Born process does not
  // exclude accompanying subprocesses with light par
  // tons).
  bool mFSet = false;
  if ( NHeavyJetVec().size() != 0 ) {
    mFSet = true;
  }

  // Partons in the initial state are massless in the CS
  // approach:
  // The following loop checks for at least one existing
  // combination 'n in addition for at least one massive
  // parton in the final state 'n for only massless par-
  // tons in the initial state.
  bool first = false;
  bool second = false;
  bool finalmass = false;
  bool initialmass = false;
  int idp = 0;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p, ++idp ) {
    if ( (*p)->coloured() && (*p)->hardProcessMass()!=ZERO && idp > 1 ) {
      finalmass = true;
    }
    if ( (*p)->coloured() && (*p)->hardProcessMass()!=ZERO && idp < 2 ) {
      initialmass = true;
    }
    if ( !first ) {
      if ( applyNotMassless(*p) ) {
        first = true;
      }
    } else {
      if ( applyNotMassless(*p) ) {
        second = true;
      }
    }
  }
  if ( first && second && (finalmass || mFSet) && !initialmass && 
       (factory()->alphaParameter() < 1.) ) {
    Repository::clog() << "DipoleMPKOperator: Warning: The alpha parameter will be set to 1.";
    Repository::clog() << "                            The massive PK Operator does not support alpha.";
    factory()->setAlphaParameter(1.);
  }

  return first && second && (finalmass || mFSet) && !initialmass;

}

bool DipoleMPKOperator::apply(tcPDPtr pd) const {
  return
    pd->hardProcessMass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

bool DipoleMPKOperator::applyNotMassless(tcPDPtr pd) const {
  return
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

void DipoleMPKOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLightJetVec().size();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLightJetVec().size();
  }
}

vector<int> DipoleMPKOperator::NHeavyJetVec() const {

  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipolePKOperator::NHeavyJetVec(): Could not find a jet particle group named 'j'" << Exception::runerror;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNHeavyJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).hardProcessMass() != ZERO )
      theNHeavyJetVec.push_back( (**theP).id() );
  }

  return theNHeavyJetVec;

}

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::me2() const {

  if ( isDR() )
    throw Exception() << "DipoleMPKOperator not implemented for dimensional reduction." << Exception::runerror;

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

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::sumParton(int id) const {

  // sumParton(int id) sums for a specific id=a'=0,1
  // over all a (see CS paper)

  // Note: We are using a parametrization where we keep s_{ja'}=2p_jp_a' fixed,
  // rather than s_{ja}=2p_jp_a, due to the substitution \eta->x/z and the sub-
  // sequent shift of the z-dependence from the hard Born ME into the PDF. 
  // Thus we need to make sure to keep the right kinematic variable fixed while
  // performig the z-integration, i.e. s_{ja'} in our case. This is partly des-
  // cribed in appendix B of the massive CS paper.
  // This also means that in the sum over heavy quark flavours, in the g->QQbar
  // contributions, we need to sum over N_F and not just N_F^{ja} (see appendix
  // B in the massive CS paper, last sentence on p. 51 in arXiv:hep-ph/0201036)
  // which is not only important for the massive PK operator here but actually 
  // also for the I operator in the massive case.

  // lambda tests the ln(m_F^2/s_ja) term between the I operator
  // (for j=gluon and k=initial state spectator) and K operator
  // (for j=gluon). If chosen coherently the physical cross
  // section has to be independent of lambda.
  double lambda=1.0;

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
    generate( ( piecewise(),
                flat(0.0,x),
                match(inverse(0.0,x,1.0) + inverse(1.0+eps,x,1.0)) ),
              r );

  z = zw.first;
  double mapz = zw.second;

  // For every new momentum fraction at which we want to evaluate a pdf
  // we introduce a new member to the pdf cache: Initialize the cache.
  // Remember, that there are potentially as many different values for 
  // z_+ (zbar_+ in the method of appendix B) as there are heavy quark
  // flavours in the jet.
  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);

  for ( map<pair<tcPDFPtr,tcPDPtr>,vector<double> >::iterator cache = 
    pdfCache.begin(); cache != pdfCache.end(); ++cache ) cache->second = nullPDFCacheVector;

  assert(pdf);

  double res = 0.;

  ////////////////////////////
  // K operator             //
  // non-color correlated   //
  ////////////////////////////

  // The non-color correlated contributions of the K operator
  // are equal in the massive and massless case.

  // Master apply prevents the DipolePKOperator though
  // from applying in the case of at least one massive
  // parton in the overall process. So need to add the 
  // corresponding contributions here (DipoleMPKOpera-
  // tor) as well.

  if ( mePartonData()[id]->id() == ParticleID::g )
    res += (KBargg() + KBarqg())*lastBorn()->me2();

  if ( abs(mePartonData()[id]->id()) < 7 )
    res += (KBarqq() + KBargq())*lastBorn()->me2();

  ////////////////////////////

  double theGammaSoft = 0.0;

  double thePqq = 0.0;
  double thePqg = 0.0;
  double thePgq = 0.0;
  double thePgg = 0.0;

  double thePqqreg = 0.0;
  double thePqgreg = 0.0;
  double thePgqreg = 0.0;
  double thePggreg = 0.0;

  double ifCorrelated = 0.0;
  double fiCorrelated = 0.0;

  int idi = 2;
  vector<Lorentz5Momentum>::const_iterator Pi = meMomenta().begin() + 2;

  /*
  double disFinite = 0.0;
  double glueFinite = 0.0;
  */

  //////////////////////////////////////////////////////////
  // Initial-Final and Final-Initial contributions        //
  //////////////////////////////////////////////////////////

  for ( cPDVector::const_iterator i = mePartonData().begin() + 2;
	i != mePartonData().end(); ++i, ++Pi, ++idi ) {

    if ( !applyNotMassless(*i) || lastBorn()->noDipole(idi,id) )
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

    //////////////////////////
    // K operator, the rest //
    //////////////////////////

    // For m_j=0 and {m_F}=empty we can use the exact massless case.
    // if ( (**i).hardProcessMass() == ZERO && lastBorn()->nHeavyJetVec().size() == 0 ) {
    if ( (**i).hardProcessMass() == ZERO && lastBorn()->nHeavyJetVec().size() == 0 ) {
      // Last term in massless K operator in (C.31) in massless paper.
      res +=
        ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
        theGammaSoft * fiCorrelated;
        // gammaSoft() * fiCorrelated;
    }

    // For m_j!=0 or {m_F}!=empty we use the massive formulae.
    else {

      Energy2 mj2 = sqr((**i).hardProcessMass());
      // Energy2 sja = 2.*( (*Pi) * meMomenta()[id]/z ); // sja=2*p_j*p_a=2*p_j*p_a'/z
      Energy2 sjaprime = 2.*( (*Pi) * meMomenta()[id] ); // sja'=2*p_j*p_a', we keep p_a' fixed during z integration

      Energy2 Qja2 = mj2 - sjaprime; // = mj2-sja' = mj2-2pjpa', which for us is constant during z integration

      // Next-to-Last term in the massive K operator in (6.55) in the massive paper.
      // Corresponds in the massless limit to the last term in the massless K operator.
      if ( mePartonData()[id]->id() == ParticleID::g ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          // ( (**i).id() == ParticleID::g ? ( Kscriptbarqg_g() + Kscriptbargg_g(Qja2) ) 
          ( (**i).id() == ParticleID::g ? ( Kscriptbarqg_g() + Kscriptbargg_g(Qja2,lambda) ) 
                                        : ( Kscriptbarqg_q(Qja2,mj2) + Kscriptbargg_q(Qja2,mj2) ) );
      }
      if ( abs(mePartonData()[id]->id()) < 7 ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          // ( (**i).id() == ParticleID::g ? ( Kscriptbarqq_g(Qja2) + Kscriptbargq_g() ) 
          ( (**i).id() == ParticleID::g ? ( Kscriptbarqq_g(Qja2,lambda) + Kscriptbargq_g() ) 
                                        : ( Kscriptbarqq_q(Qja2,mj2) + Kscriptbargq_q() ) );
      }

      // The regular splitting functions, not
      // folded with 1/z*PDF(x/z)*\Theta(z-x)
      if (thePqqreg == 0.0) thePqqreg = -1.*CF*(1.+z);
      if (thePqgreg == 0.0) thePqgreg = CF*(1.+(1.-z)*(1.-z))/z;
      if (thePgqreg == 0.0) thePgqreg = 1./2.*(z*z + (1.-z)*(1.-z));
      if (thePggreg == 0.0) thePggreg = 2.*CA*((1.-z)/z - 1. + z*(1.-z));
      // double thePqqreg = -1.*CF*(1.+z);
      // double thePqgreg = CF*(1.+(1.-z)*(1.-z))/z;
      // double thePgqreg = 1./2.*(z*z + (1.-z)*(1.-z));
      // double thePggreg = 2.*CA*((1.-z)/z - 1. + z*(1.-z));

      // Last term in massive K operator in (6.55) in massive paper.
      // Vanishes theoretically in the massless limit.
      if ( mePartonData()[id]->id() == ParticleID::g ) {
        double quarkpdfsum = 0.0;
        for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
          quarkpdfsum += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]));
          quarkpdfsum += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]));
        }
        res -=
          ifCorrelated * 
          ( ( z>x ? 1./z*( thePqgreg*quarkpdfsum + thePggreg*PDFxByz(parton) )*log((1.-z)*(sjaprime/z)/((1.-z)*(sjaprime/z)+mj2)) : 0. )
          + (gammaGluon)*PDFx(parton)*( log((sjaprime-2.*sqrt(mj2)*sqrt(sjaprime+mj2)+2.*mj2)/sjaprime) + 2.*sqrt(mj2)/(sqrt(sjaprime+mj2)+sqrt(mj2)) ) 
          );
      }
      else if ( abs(mePartonData()[id]->id()) < 7 ) {
        res -=
          ifCorrelated * 
          ( ( z>x ? 1./z*( thePqqreg*PDFxByz(parton) + thePgqreg*PDFxByz(getParticleData(ParticleID::g)) )*log((1.-z)*(sjaprime/z)/((1.-z)*(sjaprime/z)+mj2)) : 0. ) 
          + (gammaQuark)*PDFx(parton)*( log((sjaprime-2.*sqrt(mj2)*sqrt(sjaprime+mj2)+2.*mj2)/sjaprime) + 2.*sqrt(mj2)/(sqrt(sjaprime+mj2)+sqrt(mj2)) ) 
          );
      }

    }

    ////////////////
    // P operator //
    ////////////////

    // The contributions of the P operator are equal
    // in the massive and massless case.

    // Master apply prevents the DipolePKOperator though 
    // from applying in the case of at least one massive
    // parton in the overall process. So need to add the 
    // corresponding contributions here (DipoleMPKOpera-
    // tor) as well.

    double theLog = log(scale/(2.*((*Pi)*meMomenta()[id])));

    // Note: In the CS paper theLog is given by 
    // \log(\mu_F^2/(2zp_Ip_a)) = \log(\mu_F^2/(2p_Ip_a'))
    // Note: The AP splitting kernels P^{aa'} contain plus
    // distributions
    // Note however: In our implementation p_Ip_a' is kept
    // fixed, so we don't have a z dependence there

    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res +=
        ( thePgg + thePqg ) * theLog * ifCorrelated;
        // ( Pgg() + Pqg() ) * theLog * ifCorrelated;
    }

    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res +=
        ( thePqq + thePgq ) * theLog * ifCorrelated;
        // ( Pqq() + Pgq() ) * theLog * ifCorrelated;
      // if ( disFinite == 0.0 && z > x ) { 
      // extra terms only present in massless dipoles, idi is spectator
      /*
      if ( disFinite == 0.0 && z > x && mePartonData()[idi]->hardProcessMass() == ZERO ) {
        disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      */
      // if ( z > x ) 
      // extra terms only present in massless dipoles, idi is spectator
      /*
      if ( z > x && mePartonData()[idi]->hardProcessMass() == ZERO )
        res -= disFinite*ifCorrelated;
      */
    }

    if ( abs(mePartonData()[idi]->id()) < 7 ) {
      // if ( disFinite == 0.0 && z > x ) { 
      // extra terms only present in massless dipoles, idi is emitter
      /*
      if ( disFinite == 0.0 && z > x && mePartonData()[idi]->hardProcessMass() == ZERO ) {
        disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      */
      // if ( z > x ) 
      // extra terms only present in massless dipoles, idi is emitter
      /*
      if ( z > x && mePartonData()[idi]->hardProcessMass() == ZERO )
        res -= disFinite*fiCorrelated;
      */
    }

    /*
    if ( mePartonData()[idi]->id() == ParticleID::g ) {
      if ( glueFinite == 0.0 && z > x ) {
        glueFinite = 2.*CA*PDFxByz(parton)*(1.+z/6.)/z;
      }
      if ( z > x )
        res -= glueFinite*fiCorrelated;
    }
    */

  } // end loop over i

  //////////////////////////////////////////////////////////
  // Initial-Initial contributions                        //
  //////////////////////////////////////////////////////////

  if ( mePartonData()[ id == 0 ? 1 : 0 ]->coloured() &&
       !lastBorn()->noDipole(id == 0 ? 0 : 1,
                             id == 0 ? 1 : 0) ) {

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

    double theLog = log(scale/(2.*(meMomenta()[0]*meMomenta()[1])));

    // Note: In the CS paper theLog is given by 
    // \log(\mu_F^2/(2zp_Ip_a)) = \log(\mu_F^2/(2p_Ip_a'))
    // Note: The AP splitting kernels P^{aa'} contain plus
    // distributions
    // Note however: In our implementation p_Ip_a' is kept
    // fixed.

    pair<int,int> corr = id == 0 ? make_pair(0,1) : make_pair(1,0);
    double iiCorrelated = lastBorn()->colourCorrelatedME2(corr);

    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res +=
        ( thePgg + thePqg ) * theLog * iiCorrelated;
        // ( Pgg() + Pqg() ) * theLog * iiCorrelated;
      res -=
        ( KTildegg() + KTildeqg() ) * iiCorrelated;
    }    

    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res +=
        ( thePqq + thePgq ) * theLog * iiCorrelated;
        // ( Pqq() + Pgq() ) * theLog * iiCorrelated;
      res -=
        ( KTildeqq() + KTildegq() ) * iiCorrelated;
    }

  }

  //////////////////////////////////////////////////////////

  return res * mapz;

}

//////////////////////////////////////////////////////////////////////

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
  for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
    res += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]));
    res += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]));
  }
  return res*factor;
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
  for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
    res += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]));
    res += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]));
  }
  return res*factor;
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
    ( CA*( sqr(pi) - 50./9. ) + (8./9.)*lastBorn()->nLightJetVec().size() ) * PDFx(parton);
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
    ( (11./6.) * CA - (1./3.)*lastBorn()->nLightJetVec().size() + 2.*CA*log(1.-x) ) * PDFx(parton);
  if ( z > x ) {
    res += 2. * CA * ( PDFxByz(parton) - z*PDFx(parton) ) / (z*(1.-z));
    res += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::Ja_gQplus(double muQ2) const {
  double res = 
    -1. * PDFx(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(1.+muQ2/z))) );
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(2.-z+muQ2))) );
  }
  return res;
}

double DipoleMPKOperator::gammaSoft2(double muQ2) const {
  double res = 
    (1./(1.-z)) * ( -1. * PDFx(parton) * log( 1./(1.+muQ2/z) ) );
  if ( z > x ) {
    res += (1./(1.-z)) * ( 1./z * PDFxByz(parton) * log( (2.-z)/(2.-z+muQ2) ) );
  }
  return res;
}

//////////////////////////////

double DipoleMPKOperator::Kscriptqq_q(Energy2 sja, Energy2 mj2) const {
  assert(abs(parton->id()) < 7);
  double muQ2 = mj2/sja;
  double res = 
    2. * softLog(parton) + 
    Ja_gQplus(muQ2) + 
    2. * gammaSoft2(muQ2) - 
    gammaQuark/CF * PDFx(parton) + 
    ( ( muQ2==0.0 ? 0.0 : muQ2/z*log(muQ2/(z+muQ2)) ) + 1./2.*(muQ2/(z+muQ2)) ) * PDFx(parton);
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) 
         * ( -2. * log(2.-z)/(1.-z) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqg_q(Energy2 sja, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double muQ2 = mj2/sja;
  double res = 0.0;
  double factor = 1./z * ( 2.*CF/CA*( muQ2==0.0 ? 0.0 : muQ2/z*log(muQ2/(1.-z+muQ2)) ) );
  for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
    res += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]));
    res += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]));
  }
  return res*factor;
}

double DipoleMPKOperator::Kscriptgq_q() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptgg_q(Energy2 sja, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);

  // The pdf folding has to be with the gluon pdf here, not with
  // the quark pdf's. Hence the contributions of Kscriptqg_q and
  // Kscriptqq_q are implemented here separately again since the
  // safeguard assert functions would trigger if we used the im-
  // plementations above. The parton is nevertheless a gluon, so
  // we can simply use PDFx...(parton) here.

  double muQ2 = mj2/sja;
  double res = 0.0;

  // CA/CF*Kscriptqg_q contribution
  double factor = 0.0;
  if ( z > x ) {
    factor += 1./z * PDFxByz(parton) * (
      2.*CF/CA*( muQ2==0.0 ? 0.0 : muQ2/z*log(muQ2/(1.-z+muQ2)) ) );
  }
  res += CA/CF * factor;

  // Kscriptqq_q contribution
  res += 2. * softLog(parton) + 
    Ja_gQplus(muQ2) + 
    2. * gammaSoft2(muQ2) - 
    gammaQuark/CF * PDFx(parton) + 
    ( ( muQ2==0.0 ? 0.0 : muQ2/z*log(muQ2/(z+muQ2)) ) + 1./2.*(muQ2/(z+muQ2)) ) * PDFx(parton);
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) 
         * ( -2. * log(2.-z)/(1.-z) );
  }

  return res;

}

//////////////////////////////

double DipoleMPKOperator::Kscriptbarqq_q(Energy2 Qja2, Energy2 mj2) const {
  assert(abs(parton->id()) < 7);

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again
  double res = Kscriptqq_q(sjamod, mj2);

  // \deltaqq*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution
  double mubarQ2 = mj2/(-Qja2);
  res += PDFx(parton)*( sqr(pi)/3. - 2.*gsl_sf_dilog(1./(1.+mubarQ2)) + 2.*gsl_sf_dilog(-mubarQ2/(1.+mubarQ2)) 
                      + 0.5*log((1.+mubarQ2)/(1.+2.*mubarQ2)) 
                      + ( mubarQ2==0.0 ? 0.0 : 0.5*mubarQ2*(2.+mubarQ2)*log((1.+mubarQ2)/mubarQ2) )
                      - mubarQ2*(1.+mubarQ2)/(1.+2.*mubarQ2) );

  return res;

}

double DipoleMPKOperator::Kscriptbarqg_q(Energy2 Qja2, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again
  double res = Kscriptqg_q(sjamod, mj2);

  // \deltaqg*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution is zero

  return res;

}

double DipoleMPKOperator::Kscriptbargq_q() const {
  assert(abs(parton->id()) < 7);
  // Kscriptgq_q contribution is zero
  // \deltagq*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution is zero
  double res = 0.0;
  return res;
}


double DipoleMPKOperator::Kscriptbargg_q(Energy2 Qja2, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again
  double res = Kscriptgg_q(sjamod, mj2);

  // \deltagg*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution
  double mubarQ2 = mj2/(-Qja2);
  res += PDFx(parton)*( sqr(pi)/3. - 2.*gsl_sf_dilog(1./(1.+mubarQ2)) + 2.*gsl_sf_dilog(-mubarQ2/(1.+mubarQ2)) 
                      + 0.5*log((1.+mubarQ2)/(1.+2.*mubarQ2)) 
                      + ( mubarQ2==0.0 ? 0.0 : 0.5*mubarQ2*(2.+mubarQ2)*log((1.+mubarQ2)/mubarQ2) )
                      - mubarQ2*(1.+mubarQ2)/(1.+2.*mubarQ2) );

  return res;

}

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::JaNS_QQ(double muQ2) const {
  double res = 
    // 10./9. * ( 1. - sqrt(1.-4.*muQ2) ) - 
    // 8./9. * muQ2 * sqrt(1.-4.*muQ2) + 
    // 4./3. * log( (1.+sqrt(1.-4.*muQ2))/2. );
    10./9. * ( 1. - ( (1.-4.*muQ2)<=0.0 ? 0.0 : sqrt(1.-4.*muQ2) ) ) - 
    8./9. * muQ2 * ( (1.-4.*muQ2)<=0.0 ? 0.0 : sqrt(1.-4.*muQ2) ) + 
    4./3. * log( ( 1. + ( (1.-4.*muQ2)<=0.0 ? 0.0 : sqrt(1.-4.*muQ2) ) ) / 2. );
  return res;
}

double DipoleMPKOperator::Ja_QQzplus(double muQ2, int F, double zplus) const {
  double res = 0.0;
  if ( z > x && z < zplus ) {
    res += ( 1./z*PDFxByz(parton) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( zplus > x && z < zplus ) {
    res -= ( 1./zplus*PDFxByzplus(parton,F,zplus) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  return res;
}

//////////////////////////////

// double DipoleMPKOperator::Kscriptqq_g(Energy2 sja) const {
double DipoleMPKOperator::Kscriptqq_g(Energy2 sja, double lambda) const {
  assert(abs(parton->id()) < 7);

  double res = -1.*gammaGluon/CA*gammaSoft();

  Energy2 sjaprime = z*sja;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->hardProcessMass() );
    double muQ2 = mF2/sja;
    double muQ2prime = mF2/sjaprime;
    double zplus = 1./(1.+4.*muQ2prime);
    // sum only over quarks which meet special condition
    // but not if method of appendix B in massive CS pa-
    // per is used (see note at the end of appendix B)
    res += 1./(2.*CA) * (
      // PDFx(parton)*( 2./3.*(log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) ) -
      PDFx(parton)*( 2./3.*(lambda*log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) ) -
      ( x<zplus ? ( 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(zplus*muQ2prime)+5./3.) - JaNS_QQ(zplus*muQ2prime) ) ) 
                : 0. ) + 
      Ja_QQzplus(muQ2,f,zplus) + 
      // PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) 
      PDFx(parton)*( 2./3. * ( (1.-4.*muQ2prime)<=0.0 ? 0.0 : sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) ) 
    );
  }

  return res;

}

double DipoleMPKOperator::Kscriptqg_g() const {
  assert(parton->id() == ParticleID::g);
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptgq_g() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}

// double DipoleMPKOperator::Kscriptgg_g(Energy2 sja) const {
double DipoleMPKOperator::Kscriptgg_g(Energy2 sja, double lambda) const {
  assert(parton->id() == ParticleID::g);

  double res = -1.*gammaGluon/CA*gammaSoft();

  Energy2 sjaprime = z*sja;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->hardProcessMass() );
    double muQ2 = mF2/sja;
    double muQ2prime = mF2/sjaprime;
    double zplus = 1./(1.+4.*muQ2prime);
    // sum only over quarks which meet special condition
    // but not if method of appendix B in massive CS pa-
    // per is used (see note at the end of appendix B)
    res += 1./(2.*CA) * (
      // PDFx(parton)*( 2./3.*(log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) ) -
      PDFx(parton)*( 2./3.*(lambda*log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) ) -
      ( x<zplus ? ( 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(zplus*muQ2prime)+5./3.) - JaNS_QQ(zplus*muQ2prime) ) ) 
                : 0. ) + 
      Ja_QQzplus(muQ2,f,zplus) + 
      // PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) 
      PDFx(parton)*( 2./3. * ( (1.-4.*muQ2prime)<=0.0 ? 0.0 : sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) ) 
    );
  }

  return res;

}

//////////////////////////////

// double DipoleMPKOperator::Kscriptbarqq_g(Energy2 Qja2) const {
double DipoleMPKOperator::Kscriptbarqq_g(Energy2 Qja2, double lambda) const {
  assert(abs(parton->id()) < 7);

  Energy2 sjamod = -Qja2/z; // Since Qja2=-z*sja (mj2=0) this gives sja again
  // double res = Kscriptqq_g(sjamod);
  double res = Kscriptqq_g(sjamod,lambda);

  // \deltaqq*\delta(zplusbar-z)*T_R/C_A*\sum_{N_F}DeltaJbaraNS_QQbar(mF2/-Qja2) contribution
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->hardProcessMass() );
    double mubarQ2 = mF2/(-Qja2); // = mF2/(z*sja)
    double zplus = 1./(1.+4.*mubarQ2); // actually zbar_plus, see appendix B in massive CS paper
  if(x<zplus) {
    res += 1./(2.*CA)*PDFxByzplus(parton,f,zplus)*1./zplus
                     *( 8./3.*mubarQ2 - 10./9. + (10./9.+16./3.*mubarQ2)/sqrt((1.+4.*mubarQ2)*(1.+4.*mubarQ2)*(1.+4.*mubarQ2))
                      + 4./3.*( (1.-2.*mubarQ2)*sqrt(1.+4.*mubarQ2) - 1. )*log( (sqrt(1.+4.*mubarQ2)+1.)/(2.*sqrt(mubarQ2)) ) );
    }
  }

  return res;

}

double DipoleMPKOperator::Kscriptbarqg_g() const {
  assert(parton->id() == ParticleID::g);
  // Kscriptqg_g contribution is zero
  // \deltaqg*\delta(zplus-z)*... contribution is zero
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptbargq_g() const {
  assert(abs(parton->id()) < 7);
  // Kscriptgq_g contribution is zero
  // \deltagq*\delta(zplus-z)*... contribution is zero
  double res = 0.0;
  return res;
}


// double DipoleMPKOperator::Kscriptbargg_g(Energy2 Qja2) const {
double DipoleMPKOperator::Kscriptbargg_g(Energy2 Qja2, double lambda) const {
  assert(parton->id() == ParticleID::g);

  Energy2 sjamod = -Qja2/z; // Since Qja2=-z*sja (mj2=0) this gives sja again
  // double res = Kscriptgg_g(sjamod);
  double res = Kscriptgg_g(sjamod,lambda);

  // \deltagg*\delta(zplusbar-z)*T_R/C_A*\sum_{N_F}DeltaJbaraNS_QQbar(mF2/-Qja2) contribution
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->hardProcessMass() );
    double mubarQ2 = mF2/(-Qja2);
    double zplus = 1./(1.+4.*mubarQ2); // actually zbar_plus, see appendix B in massive CS paper
    if (x<zplus) {
      res += 1./(2.*CA)*PDFxByzplus(parton,f,zplus)*1./zplus
                       *( 8./3.*mubarQ2 - 10./9. + (10./9.+16./3.*mubarQ2)/sqrt((1.+4.*mubarQ2)*(1.+4.*mubarQ2)*(1.+4.*mubarQ2))
                        + 4./3.*( (1.-2.*mubarQ2)*sqrt(1.+4.*mubarQ2) - 1. )*log( (sqrt(1.+4.*mubarQ2)+1.)/(2.*sqrt(mubarQ2)) ) );
    }
  }

  return res;

}

//////////////////////////////////////////////////////////////////////

// For every new momentum fraction at which we want to evaluate a pdf
// we introduce a new member to the pdf cache.

double DipoleMPKOperator::PDFx(tcPDPtr pd) const {

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);

  map<pair<tcPDFPtr,tcPDPtr>,vector<double> >::iterator cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = nullPDFCacheVector;
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  // The convention of the pdf sets is always to return a*f(a). Upon usage, 
  // we have to remember that and rescale the result of the pdf accordingly
  // again, i.e. a*f(a)/a.
  if ( cached->second.at(0) == 0.0 ) cached->second.at(0) = pdf->xfx(particle,pd,scale,x)/x;
  return cached->second.at(0);

}

double DipoleMPKOperator::PDFxByz(tcPDPtr pd) const {

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);

  map<pair<tcPDFPtr,tcPDPtr>,vector<double> >::iterator cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = nullPDFCacheVector;
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  // The convention of the pdf sets is always to return a*f(a). Upon usage, 
  // we have to remember that and rescale the result of the pdf accordingly
  // again, i.e. a*f(a)/a.
  if ( cached->second.at(1) == 0.0 ) cached->second.at(1) = pdf->xfx(particle,pd,scale,x/z)*z/x;
  return cached->second.at(1);

}

//////////////////////////////

// For every new momentum fraction at which we want to evaluate a pdf
// we introduce a new member to the pdf cache.
// Remember, that there are potentially as many different values for 
// z_+ as there are heavy quark flavours in the jet.

double DipoleMPKOperator::PDFxByzplus(tcPDPtr pd, int F, double zplus) const {

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);

  int pdfFCacheID = 1+F;

  map<pair<tcPDFPtr,tcPDPtr>,vector<double> >::iterator cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = nullPDFCacheVector;
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  // The convention of the pdf sets is always to return a*f(a). Upon usage, 
  // we have to remember that and rescale the result of the pdf accordingly
  // again, i.e. a*f(a)/a.
  if ( cached->second.at(pdfFCacheID) == 0.0 ) cached->second.at(pdfFCacheID) = pdf->xfx(particle,pd,scale,x/zplus)*zplus/x;
  return cached->second.at(pdfFCacheID);

}

//////////////////////////////////////////////////////////////////////

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

  DipoleRepository::registerInsertionPKOperator<0,DipoleMPKOperator>("MassivePKOperator");

}

