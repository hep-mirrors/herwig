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

#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

#include <gsl/gsl_sf_dilog.h>

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
    if ( (*p)->coloured() && (*p)->mass()!=ZERO && idp > 1 ) {
      finalmass = true;
    }
    if ( (*p)->coloured() && (*p)->mass()!=ZERO && idp < 2 ) {
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

  return first && second && (finalmass || mFSet) && !initialmass;

}

bool DipoleMPKOperator::apply(tcPDPtr pd) const {
  return
    pd->mass() == ZERO &&
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
    gammaGluon = (11./6.)*CA - (1./3.)*NLightJetVec().size();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*NLightJetVec().size();
  }
}

//////////////////////////////////////////////////////////////////////

vector<int> DipoleMPKOperator::NLightJetVec() const {

  const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipolePKOperator::NLightJetVec(): Could not find a jet particle group named 'j'" << Exception::abortnow;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNLightJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).mass() == ZERO )
      theNLightJetVec.push_back( (**theP).id() );
  }

  return theNLightJetVec;

}

vector<int> DipoleMPKOperator::NHeavyJetVec() const {

  const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipolePKOperator::NHeavyJetVec(): Could not find a jet particle group named 'j'" << Exception::abortnow;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNHeavyJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).mass() != ZERO )
      theNHeavyJetVec.push_back( (**theP).id() );
  }

  return theNHeavyJetVec;

}

vector<int> DipoleMPKOperator::NLightBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNLightBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    if ( abs((**j).id()) < 7 && (**j).mass() == ZERO )
      theNLightBornVec.push_back( (**j).id() );
  }

  return theNLightBornVec;

}

vector<int> DipoleMPKOperator::NHeavyBornVec() const {

  // For the moment just count all quark and antiquark
  // constituents in the Born process.

  vector<int> theNHeavyBornVec;

  for ( cPDVector::const_iterator j = mePartonData().begin();
	j != mePartonData().end(); ++j ) {
    if ( abs((**j).id()) < 7 && (**j).mass() != ZERO )
      theNHeavyBornVec.push_back( (**j).id() );
  }

  return theNHeavyBornVec;

}

vector<int> DipoleMPKOperator::NLightProtonVec() const {

  const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("p");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipoleMPKOperator::NLightProtonVec(): Could not find a proton particle group named 'p'" << Exception::abortnow;

  const PDVector& theProtonConstitutents = theIt->second;
  vector<int> theNLightProtonVec;

  for ( PDVector::const_iterator theP = theProtonConstitutents.begin();
        theP != theProtonConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).mass() == ZERO )
      theNLightProtonVec.push_back( (**theP).id() );
  }

  return theNLightProtonVec;

}

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::me2() const {

  if ( isDR() )
    throw InitException() << "DipoleMPKOperator not implemented for dimensional reduction.";

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

  bool appendixB = true;
  // Note: We are using a parametrization where we keep s_{ja'}=2p_jp_a' fixed,
  // rather than s_{ja}=2p_jp_a, due to the substitution \eta->x/z and the sub-
  // sequent shift of the z-dependence from the hard Born ME into the PDF. 
  // Thus we need to make sure to keep the right kinematic variable fixed while
  // performig the z-integration, i.e. s_{ja'} in our case. This is partly des-
  // cribed in appendix B of the massive CS paper, but also in the last term of
  // eq. (6.55) in the massive CS paper we need to consider that s_{ja'} is our
  // fixed variable and not s_{ja}.

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

  // For every new momentum fraction at which we want to evaluate a pdf
  // we introduce a new member to the pdf cache: Initialize the cache.
  // Remember, that there are potentially as many different values for 
  // z_+ (zbar_+ in the method of appendix B) as there are heavy quark
  // flavours in the jet.
  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
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

  double disFinite = 0.0;
  double glueFinite = 0.0;

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
    // This case, however, should actually not occur here.
    if ( (**i).mass() == ZERO && NHeavyJetVec().size() == 0 ) {
      throw InitException() << "DipoleMPKOperator::sumParton: m_j=0 and {m_F}=empty. This should not have happened.";
      // Last term in massless K operator in (C.31) in massless paper.
      res +=
        ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
        theGammaSoft * fiCorrelated;
        // gammaSoft() * fiCorrelated;
    }

    // For m_j!=0 or {m_F}!=empty we use the massive formulae.
    else {

      Energy2 mj2 = sqr((**i).mass());
      Energy2 sja = 2.*( (*Pi) * meMomenta()[id]/z );

      Energy2 Qja2 = mj2 - z*sja;

      // Next-to-Last term in massive K operator in (6.55) in massive paper.
      // Corresponds in massless limit to the last term in massless K operator.
      if ( mePartonData()[id]->id() == ParticleID::g ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          ( appendixB 
          ? ( (**i).id() == ParticleID::g ? ( Kscriptbarqg_g() + Kscriptbargg_g(Qja2,appendixB) ) 
	                                  : ( Kscriptbarqg_q(Qja2,mj2) + Kscriptbargg_q(Qja2,mj2,appendixB) ) )
          : ( (**i).id() == ParticleID::g ? ( Kscriptqg_g() + Kscriptgg_g(sja,appendixB) ) 
	                                  : ( Kscriptqg_q(sja,mj2) + Kscriptgg_q(sja,mj2,appendixB) ) ) 
          );
      }
      if ( abs(mePartonData()[id]->id()) < 7 ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          ( appendixB 
          ? ( (**i).id() == ParticleID::g ? ( Kscriptbarqq_g(Qja2,appendixB) + Kscriptbargq_g() ) 
	                                  : ( Kscriptbarqq_q(Qja2,mj2,appendixB) + Kscriptbargq_q() ) )
          : ( (**i).id() == ParticleID::g ? ( Kscriptqq_g(sja,appendixB) + Kscriptgq_g() ) 
	                                  : ( Kscriptqq_q(sja,mj2,appendixB) + Kscriptgq_q() ) )
          );
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
      if (!appendixB) {
        if ( mePartonData()[id]->id() == ParticleID::g ) {
          double quarkpdfsum = 0.0;
          int nlp = NLightProtonVec().size();
          for ( int f = -nlp; f <= nlp; ++f ) {
            if ( f == 0 ) continue;
            quarkpdfsum += PDFxByz(getParticleData(f));
          }
          res -=
            ifCorrelated * 
            ( ( z>x ? 1./z*( thePqgreg*quarkpdfsum + thePggreg*PDFxByz(parton) )*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
              (gammaGluon)*PDFx(parton)*( log((sja-2.*sqrt(mj2)*sqrt(sja+mj2)+2.*mj2)/sja) + 2.*sqrt(mj2)/(sqrt(sja+mj2)+sqrt(mj2)) ) );
        }
        if ( abs(mePartonData()[id]->id()) < 7 ) {
          res -=
            ifCorrelated * 
            ( ( z>x ? 1./z*( thePqqreg*PDFxByz(parton) + thePgqreg*PDFxByz(getParticleData(ParticleID::g)) )*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
              (gammaQuark)*PDFx(parton)*( log((sja-2.*sqrt(mj2)*sqrt(sja+mj2)+2.*mj2)/sja) + 2.*sqrt(mj2)/(sqrt(sja+mj2)+sqrt(mj2)) ) );
        }
      }
      else if (appendixB) {
        if ( mePartonData()[id]->id() == ParticleID::g ) {
          double quarkpdfsum = 0.0;
          int nlp = NLightProtonVec().size();
          for ( int f = -nlp; f <= nlp; ++f ) {
            if ( f == 0 ) continue;
            quarkpdfsum += PDFxByz(getParticleData(f));
          }
          res -=
            ifCorrelated * 
            ( ( z>x ? 1./z*( thePqgreg*quarkpdfsum + thePggreg*PDFxByz(parton) )*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
              (gammaGluon)*PDFx(parton)*( log((z*sja-2.*sqrt(mj2)*sqrt(z*sja+mj2)+2.*mj2)/(z*sja)) + 2.*sqrt(mj2)/(sqrt(z*sja+mj2)+sqrt(mj2)) ) );
        }
        if ( abs(mePartonData()[id]->id()) < 7 ) {
          res -=
            ifCorrelated * 
            ( ( z>x ? 1./z*( thePqqreg*PDFxByz(parton) + thePgqreg*PDFxByz(getParticleData(ParticleID::g)) )*log((1.-z)*sja/((1.-z)*sja+mj2)) : 0. ) + 
              (gammaQuark)*PDFx(parton)*( log((z*sja-2.*sqrt(mj2)*sqrt(z*sja+mj2)+2.*mj2)/(z*sja)) + 2.*sqrt(mj2)/(sqrt(z*sja+mj2)+sqrt(mj2)) ) );
        }
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
      if ( disFinite == 0.0 && z > x && mePartonData()[idi]->mass() == ZERO ) {
        disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      // if ( z > x ) 
      // extra terms only present in massless dipoles, idi is spectator
      if ( z > x && mePartonData()[idi]->mass() == ZERO )
        res -= disFinite*ifCorrelated;
    }

    if ( abs(mePartonData()[idi]->id()) < 7 ) {
      // if ( disFinite == 0.0 && z > x ) { 
      // extra terms only present in massless dipoles, idi is emitter
      if ( disFinite == 0.0 && z > x && mePartonData()[idi]->mass() == ZERO ) {
        disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      // if ( z > x ) 
      // extra terms only present in massless dipoles, idi is emitter
      if ( z > x && mePartonData()[idi]->mass() == ZERO )
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
    // fixed, so we don't have a z dependence there

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
  int nlp = NLightProtonVec().size();
  for ( int f = -nlp; f <= nlp; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f));
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
  int nlp = NLightProtonVec().size();
  for ( int f = -nlp; f <= nlp; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f));
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
    ( CA*( sqr(pi) - 50./9. ) + (8./9.)*NLightJetVec().size() ) * PDFx(parton);
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
    ( (11./6.) * CA - (1./3.)*NLightJetVec().size() + 2.*CA*log(1.-x) ) * PDFx(parton);
  if ( z > x ) {
    res += 2. * CA * ( PDFxByz(parton) - z*PDFx(parton) ) / (z*(1.-z));
    res += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::Ja_gQplus(double muQ2) const {
  double res = 
    -1. * PDFx(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(1.+muQ2))) );
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(2.-z+muQ2))) );
  }
  return res;
}

double DipoleMPKOperator::gammaSoft2(double muQ2) const {
  double res = 
    (1./(1.-z)) * ( -1. * PDFx(parton) * log( 1./(1.+muQ2) ) );
  if ( z > x ) {
    res += (1./(1.-z)) * ( 1./z * PDFxByz(parton) * log( (2.-z)/(2.-z+muQ2) ) );
  }
  return res;
}

double DipoleMPKOperator::Ja_gQplusmod(double muQ2) const {
  double res = 
    -1. * PDFx(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(1.+muQ2/z))) );
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * ( (1.-z)/(2.*(1.-z+muQ2)*(1.-z+muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(2.-z+muQ2))) );
  }
  return res;
}

double DipoleMPKOperator::gammaSoft2mod(double muQ2) const {
  double res = 
    (1./(1.-z)) * ( -1. * PDFx(parton) * log( 1./(1.+muQ2/z) ) );
  if ( z > x ) {
    res += (1./(1.-z)) * ( 1./z * PDFxByz(parton) * log( (2.-z)/(2.-z+muQ2) ) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqq_q(Energy2 sja, Energy2 mj2, bool appendixB) const {
  assert(abs(parton->id()) < 7);
  double muQ2 = mj2/sja;
  double res = 
    2. * softLog(parton) + 
    ( appendixB ? Ja_gQplusmod(muQ2) : Ja_gQplus(muQ2) ) + 
    2. * ( appendixB ? gammaSoft2mod(muQ2) : gammaSoft2(muQ2) ) - 
    gammaQuark/CF * PDFx(parton) + 
    ( appendixB ? muQ2/z*log(muQ2/(z+muQ2)) + 1./2.*(muQ2/(z+muQ2)) : muQ2*log(muQ2/(1.+muQ2)) + 1./2.*(muQ2/(1.+muQ2)) ) * PDFx(parton);
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * (
      -2. * log(2.-z)/(1.-z) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqg_q(Energy2 sja, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double muQ2 = mj2/sja;
  double res = 0.0;
  double factor = 1./z * ( 2.*CF/CA*muQ2/z*log(muQ2/(1.-z+muQ2)) );
  int nlp = NLightProtonVec().size();
  for ( int f = -nlp; f <= nlp; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f));
  }
  return res*factor;
}

double DipoleMPKOperator::Kscriptgq_q() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptgg_q(Energy2 sja, Energy2 mj2, bool appendixB) const {
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
      2.*CF/CA*muQ2/z*log(muQ2/(1.-z+muQ2)) );
  }
  res += CA/CF * factor;

  // Kscriptqq_q contribution
  res += 2. * softLog(parton) + 
    ( appendixB ? Ja_gQplusmod(muQ2) : Ja_gQplus(muQ2) ) + 
    2. * ( appendixB ? gammaSoft2mod(muQ2) : gammaSoft2(muQ2) ) - 
    gammaQuark/CF * PDFx(parton) + 
    ( appendixB ? muQ2/z*log(muQ2/(z+muQ2)) + 1./2.*(muQ2/(z+muQ2)) : muQ2*log(muQ2/(1.+muQ2)) + 1./2.*(muQ2/(1.+muQ2)) ) * PDFx(parton);
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * (
      -2. * log(2.-z)/(1.-z) );
  }

  return res;

}

double DipoleMPKOperator::Kscriptbarqq_q(Energy2 Qja2, Energy2 mj2, bool appendixB) const {
  assert(abs(parton->id()) < 7);

  double mubarQ2 = mj2/(-Qja2);
  double res = 0.0;

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again

  res += Kscriptqq_q(sjamod, mj2, appendixB)
      // \deltaqq*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution
      + PDFx(parton)*( sqr(pi)/3. - 2.*gsl_sf_dilog(1./(1.+mubarQ2)) + 2.*gsl_sf_dilog(-mubarQ2/(1.+mubarQ2)) 
                     + 0.5*log((1.+mubarQ2)/(1.+2.*mubarQ2)) + 0.5*mubarQ2*(2.+mubarQ2)*log((1.+mubarQ2)/mubarQ2) 
                     - mubarQ2*(1.+mubarQ2)/(1.+2.*mubarQ2) );

  return res;

}

double DipoleMPKOperator::Kscriptbarqg_q(Energy2 Qja2, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);

  double res = 0.0;

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again

  res += Kscriptqg_q(sjamod, mj2);
  // \deltaqg*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution is zero

  return res;

}

double DipoleMPKOperator::Kscriptbargq_q() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}


double DipoleMPKOperator::Kscriptbargg_q(Energy2 Qja2, Energy2 mj2, bool appendixB) const {
  assert(parton->id() == ParticleID::g);

  double mubarQ2 = mj2/(-Qja2);
  double res = 0.0;

  Energy2 sjamod = (mj2-Qja2)/z; // Since Qja2=mj2-z*sja this gives sja again

  res += Kscriptgg_q(sjamod, mj2, appendixB)
      // \deltagg*\delta(1-z)*DeltaJbaraNS_gQ(mj2/-Qja2) contribution
      + PDFx(parton)*( sqr(pi)/3. - 2.*gsl_sf_dilog(1./(1.+mubarQ2)) + 2.*gsl_sf_dilog(-mubarQ2/(1.+mubarQ2)) 
                     + 0.5*log((1.+mubarQ2)/(1.+2.*mubarQ2)) + 0.5*mubarQ2*(2.+mubarQ2)*log((1.+mubarQ2)/mubarQ2) 
                     - mubarQ2*(1.+mubarQ2)/(1.+2.*mubarQ2) );

  return res;

}

//////////////////////////////

double DipoleMPKOperator::JaNS_QQ(double muQ2) const {
  double res = 
    10./9. * ( 1. - sqrt(1.-4.*muQ2) ) - 
    8./9. * muQ2 * sqrt(1.-4.*muQ2) + 
    4./3. * log( (1.+sqrt(1.-4.*muQ2))/2. );
  return res;
}

double DipoleMPKOperator::Ja_QQzplus(double muQ2, int F, double zplus) const {
  double res = 0.0;
  if ( z > x && z < zplus ) {
    res += ( 1./z*PDFxByz(parton) - 1./zplus*PDFxByzplus(parton,F,zplus) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( z > x && z < zplus ) {
    res += ( 1./zplus*PDFxByzplus(parton,F,zplus) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( zplus > x && z < zplus ) {
    res -= ( 1./zplus*PDFxByzplus(parton,F,zplus) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  return res;
}

double DipoleMPKOperator::Ja_QQzplusmod(double muQ2, int F, double zplus) const {
// This should actually be the same as Ja_QQzplus()
  double res = 0.0;
  if ( z > x && z < zplus ) {
    res += ( 1./z*PDFxByz(parton) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( zplus > x && z < zplus ) {
    res -= ( 1./zplus*PDFxByzplus(parton,F,zplus) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqq_g(Energy2 sja, bool appendixB) const {
  assert(abs(parton->id()) < 7);

  double res = -1.*gammaGluon/CA*gammaSoft();

  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(NHeavyJetVec()[f])->mass() );
    double muQ2 = mF2/sja;
    double muQ2prime = mF2/(z*sja);
    double muQ2bar = mF2/(-z*sja); // only for modified zplus, see appendix B, remember that m_j=m_g=0
    double zplus = (appendixB ? 1./(1.+4.*muQ2bar) : 1. - 4*muQ2);
    // sum only over heavy quarks which meet special condition
    // but not if method of appendix B is used (see comment at
    // the end of appendix B, also need different zplus then)
    if( !appendixB && sja <= 4.*mF2 ) continue;
    res += 1./(2.*CA) * (
      PDFx(parton)*( ( appendixB ? 2./3.*(log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) 
                                 : 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) ) -
      ( x<zplus ? ( appendixB ? 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(zplus*muQ2prime)+5./3.) - JaNS_QQ(zplus*muQ2prime) ) 
                              : 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) ) : 0. ) + 
      ( appendixB ? Ja_QQzplusmod(muQ2,f,zplus) : Ja_QQzplus(muQ2,f,zplus) ) + 
      ( appendixB ? PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) 
                  : PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) ) 
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

double DipoleMPKOperator::Kscriptgg_g(Energy2 sja, bool appendixB) const {
  assert(parton->id() == ParticleID::g);

  double res = -1.*gammaGluon/CA*gammaSoft();

  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(NHeavyJetVec()[f])->mass() );
    double muQ2 = mF2/sja;
    double muQ2prime = mF2/(z*sja);
    double muQ2bar = mF2/(-z*sja); // only for modified zplus, see appendix B, remember that m_j=m_g=0
    double zplus = (appendixB ? 1./(1.+4.*muQ2bar) : 1. - 4*muQ2);
    // sum only over heavy quarks which meet special condition
    // but not if method of appendix B is used (see comment at
    // the end of appendix B, also need different zplus then)
    if( !appendixB && sja <= 4.*mF2 ) continue;
    res += 1./(2.*CA) * (
      PDFx(parton)*( ( appendixB ? 2./3.*(log(muQ2prime)+5./3.) - JaNS_QQ(muQ2prime) 
                                 : 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) ) -
      ( x<zplus ? ( appendixB ? 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(zplus*muQ2prime)+5./3.) - JaNS_QQ(zplus*muQ2prime) ) 
                              : 1./zplus*PDFxByzplus(parton,f,zplus)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) ) : 0. ) + 
      ( appendixB ? Ja_QQzplusmod(muQ2,f,zplus) : Ja_QQzplus(muQ2,f,zplus) ) + 
      ( appendixB ? PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2prime)*(1.-4.*muQ2prime)*(1.-4.*muQ2prime)) ) 
                  : PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) ) 
      );
  }

  return res;

}

double DipoleMPKOperator::Kscriptbarqq_g(Energy2 Qja2, bool appendixB) const {
  assert(abs(parton->id()) < 7);

  double res = 0.0;

  Energy2 sjamod = -Qja2/z; // Since Qja2=-z*sja (mj2=0) this gives sja again
  res += Kscriptqq_g(sjamod,appendixB); // z_+ is automatically modified

  // + \deltaqq*\delta(zplusbar-z)*T_R/C_A*\sum_{N_F}DeltaJbaraNS_QQbar(mF2/-Qja2) contribution
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(NHeavyJetVec()[f])->mass() );
    double mubarQ2 = mF2/(-Qja2);
    double zplus = 1./(1.+4.*mubarQ2); // actually zbar_plus, see appendix B
    res += 1./(2.*CA)*PDFxByzplus(parton,f,zplus)
                     *( 8./3.*mubarQ2 - 10./9. + (10./9.+16./3.*mubarQ2)/sqrt((1.+4.*mubarQ2)*(1.+4.*mubarQ2)*(1.+4.*mubarQ2))
                      + 4./3.*( (1.-2.*mubarQ2)*sqrt(1.+4.*mubarQ2) - 1. )*log( (sqrt(1.+4.*mubarQ2)+1.)/(2.*sqrt(mubarQ2)) ) );
  }

  return res;

}

double DipoleMPKOperator::Kscriptbarqg_g() const {
  assert(parton->id() == ParticleID::g);
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptbargq_g() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}


double DipoleMPKOperator::Kscriptbargg_g(Energy2 Qja2, bool appendixB) const {
  assert(parton->id() == ParticleID::g);

  double res = 0.0;

  Energy2 sjamod = -Qja2/z; // Since Qja2=-z*sja (mj2=0) this gives sja again
  res += Kscriptgg_g(sjamod,appendixB); // z_+ is automatically modified

  // + \deltagg*\delta(zplusbar-z)*T_R/C_A*\sum_{N_F}DeltaJbaraNS_QQbar(mF2/-Qja2) contribution
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(NHeavyJetVec()[f])->mass() );
    double mubarQ2 = mF2/(-Qja2);
    double zplus = 1./(1.+4.*mubarQ2); // actually zbar_plus, see appendix B
    res += 1./(2.*CA)*PDFxByzplus(parton,f,zplus)
                     *( 8./3.*mubarQ2 - 10./9. + (10./9.+16./3.*mubarQ2)/sqrt((1.+4.*mubarQ2)*(1.+4.*mubarQ2)*(1.+4.*mubarQ2))
                      + 4./3.*( (1.-2.*mubarQ2)*sqrt(1.+4.*mubarQ2) - 1. )*log( (sqrt(1.+4.*mubarQ2)+1.)/(2.*sqrt(mubarQ2)) ) );
  }

  return res;

}

//////////////////////////////////////////////////////////////////////

// For every new momentum fraction at which we want to evaluate a pdf
// we introduce a new member to the pdf cache.

double DipoleMPKOperator::PDFx(tcPDPtr pd) const {

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
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
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
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
// z_+ (zbar_+ in the method of appendix B) as there are heavy quark
// flavours in the jet.

double DipoleMPKOperator::PDFxByzplus(tcPDPtr pd, int F, double zplus) const {

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<NHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
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

  DipoleRepository::registerInsertionOperator<0,DipoleMPKOperator>("MassivePKOperator");

}

