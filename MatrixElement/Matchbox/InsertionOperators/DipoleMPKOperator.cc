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

void DipoleMPKOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    // gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLight();
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLightJetVec().size();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    // KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLight();
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLightJetVec().size();
  }
}

bool DipoleMPKOperator::apply(tcPDPtr pd) const {
  return
    pd->mass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

bool DipoleMPKOperator::applyNotMassless(tcPDPtr pd) const {
  return
    // pd->mass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}

bool DipoleMPKOperator::apply(const cPDVector& pd) const {

  if ( !apply(pd[0]) && !apply(pd[1]) ) {
    cout << "DipoleMPKOperator::apply (master apply): ( !apply(pd[0]) && !apply(pd[1]) ). Return false!" << endl;
    return false;
  }

  // DipoleMIOperator should only apply, if massive
  // partons exist in the given process (aka in the
  // final state):
  // This applies also to gluons in the final state of
  // the LO process and further splittings g->Q\bar{Q}
  // inside the jet.
  bool mFSet = false;
  if ( NHeavyJetVec().size() != 0 ) {
    mFSet = true;
  }

  // Partons in the initial state are assumed massless in
  // the CS approach.
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

  if ( first && second && (finalmass || mFSet) && !initialmass ) {
    cout << "DipoleMPKOperator::apply (master apply): Return true!" << endl;
    cout << endl;
    cout << "     !!!!! Attention !!!!!" << endl;
    cout << "     Number of massless flavours in jet particle group (aka n_f) = " << NLightJetVec().size() << endl;
    cout << "     Number of massive flavours in jet particle group (aka n_F or n_{f,h}) = " << NHeavyJetVec().size() << endl;
    cout << "     Ensure consistent usage!" << endl;
    cout << endl;
  }
  if ( !( first && second && (finalmass || mFSet) && !initialmass ) )
    cout << "DipoleMPKOperator::apply (master apply): Return false!" << endl;

  return first && second && (finalmass || mFSet) && !initialmass;

}

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
    // if ( (**j).id() > 0 && (**j).id() < 7 && (**j).mass() == ZERO )
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
    // if ( (**j).id() > 0 && (**j).id() < 7 && (**j).mass() != ZERO )
    if ( abs((**j).id()) < 7 && (**j).mass() != ZERO )
      theNHeavyBornVec.push_back( (**j).id() );
  }

  return theNHeavyBornVec;

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

    //////////////////////////
    // K operator, the rest //
    //////////////////////////

    // For m_j=0 and {m_F}=empty we can use the exact massless case.
    // This case, however, should actually not occur here.
    if ( (**i).mass() == ZERO && lastBorn()->nHeavyJetVec().size() == 0 ) {
      throw InitException() << "DipoleMPKOperator::sumParton: m_j=0 and {m_F}=empty. This should not have happened.";
      // Last term in massless K operator in (C.31) in massless paper.
      res +=
        ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
        theGammaSoft * fiCorrelated;
    }

    // For m_j!=0 or {m_F}!=empty we use the massive formulae.
    else {

      Energy2 mj2 = sqr((**i).mass());
      Energy2 sja = 2.*( (*Pi) * meMomenta()[id]/z );

      // Next-to-Last term in massive K operator in (6.55) in massive paper.
      // Corresponds in massless limit to the last term in massless K operator.
      if ( mePartonData()[id]->id() == ParticleID::g ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          ( (**i).id() == ParticleID::g ? ( Kscriptqg_g() + Kscriptgg_g(sja) ) : ( Kscriptqg_q(sja,mj2) + Kscriptgg_q(sja,mj2) ) );
      }
      if ( abs(mePartonData()[id]->id()) < 7 ) {
        res -=
          ( (**i).id() == ParticleID::g ? CA : CF ) *
          fiCorrelated * 
          ( (**i).id() == ParticleID::g ? ( Kscriptqq_g(sja) + Kscriptgq_g() ) : ( Kscriptqq_q(sja,mj2) + Kscriptgq_q() ) );
      }

      // The regular splitting functions, not
      // folded with 1/z*PDF(x/z)*\Theta(z-x)
      double thePqqreg = -1.*CF*(1.+z);
      double thePqgreg = CF*(1.+(1.-z)*(1.-z))/z;
      double thePgqreg = 1./2.*(z*z + (1.-z)*(1.-z));
      double thePggreg = 2.*CA*((1.-z)/z - 1. + z*(1.-z));

      // Last term in massive K operator in (6.55) in massive paper.
      // Vanishes theoretically in the massless limit.
      if ( mePartonData()[id]->id() == ParticleID::g ) {
        double quarkpdfsum = 0.0;
        int nl= lastBorn()->nLight();
        for ( int f = -lastBorn()->nLight(); f <= nl; ++f ) {
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

    ////////////////
    // P operator //
    ////////////////

    // The contributions of the P operator are equal
    // in the massive and massless case.

    // The extra terms, which render the dipole kernels
    // positive definite, are subtracted here again in
    // integrated form (disfinite and gluefinite).

    // NOTE: At the moment the same extra terms are also
    // in use in the massive dipoles. This has yet to be 
    // clarified, though, for the massive case.

    double theLog = log(scale/(2.*((*Pi)*meMomenta()[id])));

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

    pair<int,int> corr = id == 0 ? make_pair(0,1) : make_pair(1,0);
    double iiCorrelated = lastBorn()->colourCorrelatedME2(corr);

    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res +=
	( thePgg + thePqg ) * theLog * iiCorrelated;
      res -=
	( KTildegg() + KTildeqg() ) * iiCorrelated;
    }    

    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res +=
	( thePqq + thePgq ) * theLog * iiCorrelated;
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

//////////////////////////////////////////////////////////////////////

double DipoleMPKOperator::Ja_gQplus(double muQ2) const {
//   assert(abs(parton->id()) < 7);
  double res = 
    -1. * PDFx(parton) * ( (1.-z)/(2.*(1.-z-muQ2)*(1.-z-muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(1.+muQ2))) );
  if ( z > x ) {
    res += 1./z * PDFxByz(parton) * ( (1.-z)/(2.*(1.-z-muQ2)*(1.-z-muQ2)) - 2./(1.-z)*(1.+log((1.-z+muQ2)/(2.-z+muQ2))) );
  }
  return res;
}

double DipoleMPKOperator::gammaSoft2(double muQ2) const {
//   assert(abs(parton->id()) < 7);
  double res = 
    (1./(1.-z)) * ( -1. * PDFx(parton) * 2. * log( 1./(1.+muQ2) ) );
  if ( z > x ) {
    res += (1./(1.-z)) * ( 1./z * PDFxByz(parton) * 2. * log( (2.-z)/(2.-z+muQ2) ) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqq_q(Energy2 sja, Energy2 mj2) const {
  assert(abs(parton->id()) < 7);
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

double DipoleMPKOperator::Kscriptqg_q(Energy2 sja, Energy2 mj2) const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double muQ2 = mj2/sja;
  double res = 0.0;
  double factor = 1./z * ( 2.*CF/CA*muQ2/z*log(muQ2/(1.-z+muQ2)) );
  int nl= lastBorn()->nLight();
  for ( int f = -lastBorn()->nLight(); f <= nl; ++f ) {
    if ( f == 0 )
      continue;
    res += PDFxByz(getParticleData(f))*factor;
  }
  return res;
}

double DipoleMPKOperator::Kscriptgq_q() const {
  assert(abs(parton->id()) < 7);
  double res = 0.0;
  return res;
}

double DipoleMPKOperator::Kscriptgg_q(Energy2 sja, Energy2 mj2) const {

  // I think the pdf folding has to be with the gluon pdf, not with the
  // quark pdf's, for both contributing terms here, because the overall
  // contribution is for (aa'_j)=(gg_q), although being a sum of the
  // (aa'_j)=(qg_q) and (aa'_j)=(qq_q) contributions. -- CR 20140403

  assert(parton->id() == ParticleID::g);
  double muQ2 = mj2/sja;
  double res = 0.0;

  // CA/CF*Kscriptqg_q(sja,mj2) contribution
  double factor = 0.0;
  if ( z > x ) {
    factor += 1./z * PDFxByz(parton) * (
      2.*CF/CA*muQ2/z*log(muQ2/(1.-z+muQ2)) );
  }
  res += CA/CF * factor;

  // Kscriptqq_q(sja,mj2) contribution
  res += 2. * softLog(parton) + 
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

//////////////////////////////

double DipoleMPKOperator::JaNS_QQ(double muQ2) const {
  double res = 
    10./9. * ( 1. - sqrt(1.-4.*muQ2) ) - 
    8./9. * muQ2 * sqrt(1.-4.*muQ2) + 
    4./3. * log( (1.+sqrt(1.-4.*muQ2))/2. );
  return res;
}

double DipoleMPKOperator::Ja_QQzplus(double muQ2, int F) const {
  double zplus = 1. - 4.*muQ2;
  double res = 0.0;
  if ( z > x && z < zplus ) {
    res += ( 1./z*PDFxByz(parton) - 1./zplus*PDFxByzplus(parton,F,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( z > x && z < zplus ) {
    res += ( 1./zplus*PDFxByzplus(parton,F,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  if ( zplus > x && z < zplus ) {
    res -= ( 1./zplus*PDFxByzplus(parton,F,muQ2) ) * 2./3.*( (1.-z+2.*muQ2)/((1.-z)*(1.-z)) * sqrt(1.-4.*muQ2/(1.-z)) );
  }
  return res;
}

double DipoleMPKOperator::Kscriptqq_g(Energy2 sja) const {
  assert(abs(parton->id()) < 7);
  double res = -1.*gammaGluon/CA*gammaSoft();
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->mass() );
    double muQ2 = mF2/sja;
    double zplus = 1. - 4*muQ2;
    if( sja <= 4.*mF2 ) continue; // sum only over heavy quarks with special condition
    res += 1./(2.*CA) * (
      PDFx(parton)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) -
      1./zplus*PDFxByzplus(parton,f,muQ2)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) + 
      Ja_QQzplus(muQ2,f) + 
      PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) );
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

double DipoleMPKOperator::Kscriptgg_g(Energy2 sja) const {
  assert(parton->id() == ParticleID::g);
  double res = -1.*gammaGluon/CA*gammaSoft();
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) { // sum over heavy flavours
    Energy2 mF2 = sqr( getParticleData(lastBorn()->nHeavyJetVec()[f])->mass() );
    double muQ2 = mF2/sja;
    double zplus = 1. - 4*muQ2;
    if( sja <= 4.*mF2 ) continue; // sum only over quarks with special condition
    res += 1./(2.*CA) * (
      PDFx(parton)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) -
      1./zplus*PDFxByzplus(parton,f,muQ2)*( 2./3.*(log(muQ2)+5./3.) - JaNS_QQ(muQ2) ) + 
      Ja_QQzplus(muQ2,f) + 
      PDFx(parton)*( 2./3.*sqrt((1.-4.*muQ2)*(1.-4.*muQ2)*(1.-4.*muQ2)) ) );
  }
  return res;
}

//////////////////////////////////////////////////////////////////////

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
  if ( cached->second.at(1) == 0.0 ) cached->second.at(1) = pdf->xfx(particle,pd,scale,x/z)*z/x;
  return cached->second.at(1);

}

//////////////////////////////

double DipoleMPKOperator::PDFxByzplus(tcPDPtr pd, int F, double muQ2) const {

  // NOTE: Might need some debugging. Not sure about the
  // multiplication by zplus/x at the end -- CR 20140425

  vector<double> nullPDFCacheVector;
  for( size_t f=0; f<lastBorn()->nHeavyJetVec().size(); ++f ) nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);
  nullPDFCacheVector.push_back(0.0);

  double zplus = 1. - 4.*muQ2;
  int pdfFCacheID = 1+F;

  map<pair<tcPDFPtr,tcPDPtr>,vector<double> >::iterator cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = nullPDFCacheVector;
    cached = pdfCache.find(make_pair(pdf,pd));
  }
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

