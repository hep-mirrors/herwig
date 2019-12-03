// -*- C++ -*-
//
// DipolePKOperator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipolePKOperator class.
//

#include "DipolePKOperator.h"
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

using namespace Herwig;
using Constants::pi;

DipolePKOperator::DipolePKOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    KQuark(-1.0), KGluon(-1.0),
    scale(ZERO), x(0.), z(0.) {}

DipolePKOperator::~DipolePKOperator() {}

IBPtr DipolePKOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipolePKOperator::fullclone() const {
  return new_ptr(*this);
}

//////////////////////////////////////////////////////////////////////

bool DipolePKOperator::apply(const cPDVector& pd) const {

  // DipolePKOperator should only apply if in the overall
  // process only massless partons can occur.

  if ( !apply(pd[0]) && !apply(pd[1]) ) {
    return false;
  }

  // Prohibit splittings g->Q\bar{Q} in the final state.
  // These are covered completely by DipoleMPKOperator.
  if ( NHeavyJetVec().size()!=0 ) {
    return false;
  }

  bool first = false;
  bool second = false;
  for ( cPDVector::const_iterator p = pd.begin();
	p != pd.end(); ++p ) {
    // Since this loop only checks for at least one exis-
    // ting combination: Return false if any massive par-
    // tons are present (covered by DipoleMPKOperator).
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

bool DipolePKOperator::apply(tcPDPtr pd) const {
  return
    pd->hardProcessMass() == ZERO &&
    (abs(pd->id()) < 7 || pd->id() == ParticleID::g);
}


void DipolePKOperator::setAlpha(double alpha)const{
  factory()->setAlphaParameter(alpha);
}


void DipolePKOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    // gammaGluon = (11./6.)*CA - (1./3.)*NLightJetVec().size();
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLightJetVec().size();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    // KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*NLightJetVec().size();
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLightJetVec().size();
    if ( isDR() ) {
      gammaQuark -= CF/2.;
      gammaGluon -= CA/6.;
    }
  }
}

//////////////////////////////////////////////////////////////////////

vector<int> DipolePKOperator::NLightJetVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("j");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipolePKOperator::NLightJetVec(): Could not find a jet particle group named 'j'" << Exception::runerror;

  const PDVector& theJetConstitutents = theIt->second;
  vector<int> theNLightJetVec;

  for ( PDVector::const_iterator theP = theJetConstitutents.begin();
        theP != theJetConstitutents.end(); ++theP ) {
    if ( (**theP).id() > 0 && (**theP).id() < 7 && (**theP).hardProcessMass() == ZERO )
      theNLightJetVec.push_back( (**theP).id() );
  }

  return theNLightJetVec;

}

vector<int> DipolePKOperator::NHeavyJetVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
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

vector<int> DipolePKOperator::NLightBornVec() const {

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

vector<int> DipolePKOperator::NHeavyBornVec() const {

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

vector<int> DipolePKOperator::NLightProtonVec() const {

  // const map<string,PDVector>& theParticleGroups = MatchboxFactory::currentFactory()->particleGroups();
  const map<string,PDVector>& theParticleGroups = factory()->particleGroups();
  map<string,PDVector>::const_iterator theIt = theParticleGroups.find("p");
  if ( theIt == theParticleGroups.end() )
    throw Exception() << "DipolePKOperator::NLightProtonVec(): Could not find a proton particle group named 'p'" << Exception::runerror;

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

double DipolePKOperator::me2() const {

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

double DipolePKOperator::sumParton(int id) const {

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

  // For every new momentum fraction at which we want to evaluate a pdf
  // we introduce a new member to the pdf cache: Initialize the cache.
  for ( map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator cache =
	  pdfCache.begin(); cache != pdfCache.end(); ++cache )
    cache->second = make_pair(0.0,0.0);

  assert(pdf);

  double res = 0.;
  

  ////////////////////////////
  // K operator             //
  // non-color correlated   //
  ////////////////////////////

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

  /*
  double disFinite = 0.0;
  double glueFinite = 0.0;
  */

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

    // Last term in massless K operator in (C.31) in massless paper
    res +=
      ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
      theGammaSoft * fiCorrelated;

    ////////////////
    // P operator //
    ////////////////

    // The extra terms, which render the dipole kernels
    // positive definite, are subtracted here again in
    // integrated form (disfinite and gluefinite).

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
    }

    if ( abs(mePartonData()[id]->id()) < 7 ) {
      res +=
	( thePqq + thePgq ) * theLog * ifCorrelated;
      /*
      if ( disFinite == 0.0 && z > x ) {
        disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      if ( z > x )
        res -= disFinite*ifCorrelated;
      */
    }

    /*
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
    // fixed, so we don't have a z dependence there

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



double DipolePKOperator::softLogByz(tcPDPtr p) const {
  double res = ( sqr(log(1.-x))/2. - sqr(pi)/6. ) * PDFx(p);
  if ( z > x ) {
    res += (PDFxByz(p) - z*PDFx(p))*log(1.-z)/(z*(1.-z));
    res -= PDFxByz(p)*log(z)/(z*(1.-z));
  }
  return res;
}

double DipolePKOperator::softLog(tcPDPtr p) const {
  double res =0.;
  if(x>(1.-factory()->alphaParameter()))
     res +=(sqr(log(1.-x)) - sqr(log(factory()->alphaParameter()))) * PDFx(p) / 2.;
  
  if(z>x&&z>(1.-factory()->alphaParameter()))
     res += (PDFxByz(p)/z  - PDFx(p))*log(1.-z)/(1.-z);
  
  return res;
}

double DipolePKOperator::gammaSoft() const {
  double res =factory()->alphaParameter()*PDFx(parton);
  if(x>(1.-factory()->alphaParameter()))
     res += ( log(1.-x)- log(factory()->alphaParameter()) )*PDFx(parton);
  if(z>x&&z>(1.-factory()->alphaParameter()))
     res += (PDFxByz(parton)/z  - PDFx(parton)) / (1.-z);
  return res;
}


double DipolePKOperator::KBarqq() const {
  assert(abs(parton->id()) < 7);
  double res =  2.* CF*softLogByz(parton);
  
  res +=  CF*(sqr(pi) - 5.)*PDFx(parton);//(C.17)
  
  if ( z > x ) {
    res += ( CF*(1.-z) - CF*(1.+z)*log(factory()->alphaParameter()*(1.-z)/z) )* PDFxByz(parton)/ z;
  }
  
  double alpha = factory()->alphaParameter();

  if ( alpha < 1. ) {
    res+=PDFx(parton)*(2.*CF*sqr(log(alpha))-gammaQuark*(alpha-1.-log(alpha)));
    if ( z > x ) {
      res+=CF*(2./(1.-z)*(log(alpha*(2.-z)/(1.+alpha-z))-log((2.-z)/(1.-z))*((z<(1.-alpha))?1.:0.)))*PDFxByz(parton)/z;
    }
  }
  
  return res ;
}

double DipolePKOperator::KBargg() const {

  assert(parton->id() == ParticleID::g);
  double res = 
    2.* CA* softLogByz(parton);
  
  res +=   ( CA*( sqr(pi) - 50./9. ) + (8./9.)*lastBorn()->nLightJetVec().size() ) * PDFx(parton);
  if ( z > x ) {
    res += 2.*CA*((1.-z)/z-1.+z*(1.-z))*log(factory()->alphaParameter()*(1.-z)/z)*PDFxByz(parton)/z;
  }
  
  if ( factory()->alphaParameter() < 1. ) {
    double alpha = factory()->alphaParameter();
    res += PDFx(parton)*(2.*CA*sqr(log(alpha))-gammaGluon*(alpha-1.-log(alpha)));
    if ( z > x ) {
      res+=CA*(2./(1.-z)*(log(alpha*(2.-z)/(1.+alpha-z))-log((2.-z)/(1.-z))*((z<(1.-alpha))?1.:0.)))*PDFxByz(parton)/z;
    }
  }
  
  
  return res;
}

double DipolePKOperator::KBargq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return
    PDFxByz(getParticleData(ParticleID::g)) *
    ( 0.5*(sqr(z)+sqr(1.-z))*log(factory()->alphaParameter()*(1.-z)/z) + z*(1.-z) ) / z;
}

double DipolePKOperator::KBarqg() const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double res = 0.0;
  double factor = CF * ( ( (1.+sqr(1.-z)) / z ) * log(factory()->alphaParameter()*(1.-z)/z) + z ) / z;
  // // int nlp = NLightProtonVec().size();
  // int nlp = lastBorn()->nLightProtonVec().size();
  // for ( int f = -nlp; f <= nlp; ++f ) {
  //   if ( f == 0 )
  //     continue;
  //   res += PDFxByz(getParticleData(f))*factor;
  // }
  for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
    res += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]))*factor;
    res += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]))*factor;
  }
  return res;
}

double DipolePKOperator::KTildeqq() const {
  assert(abs(parton->id()) < 7);
  double res =
    2.*CF*softLog(parton) - CF*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res -= ( CF * (1.+z) * log((1.-z)/factory()->alphaParameter()) ) * PDFxByz(parton) / z;
  }
  
  if ( factory()->alphaParameter() != 1. ) {
    double alpha=factory()->alphaParameter();
    if ( z > x ) {
      double AP_Pqqhat=CF *((1.+z*z)/(1.-z));
      res+=AP_Pqqhat*log(min(1.,alpha/(1.-z))) *PDFxByz(parton)/z;
      res+= (CF*2./(1.-z)*(log((1.+alpha-z)/alpha)-((z>(1.-alpha))?1.:0.)*log(2.-z))   *PDFxByz(parton)/z);
    }
  }
  
  return res;
}

double DipolePKOperator::KTildegg() const {
  
  assert(parton->id() == ParticleID::g);
  double res =
    2.*CA*softLog(parton) - CA*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res += ( 2.*CA * ( (1.-z)/z -1. + z*(1.-z) ) * log((1.-z)/factory()->alphaParameter()) ) * PDFxByz(parton) / z;
  }
  double alpha = factory()->alphaParameter();
  if ( alpha != 1. ) {
    if ( z > x ) {
      double AP_Pqqhat=2.*CA *((z)/(1.-z)+(1.-z)/(z)+z*(1.-z));
      res+=AP_Pqqhat*log(min(1.,alpha/(1.-z))) *PDFxByz(parton)/z;
      res+= (CA*2./(1.-z)*(log((1.+alpha-z)/alpha)-((z>(1.-alpha))?1.:0.)*log(2.-z))   *PDFxByz(parton)/z);
    }
  }
  
  return res;
}

double DipolePKOperator::KTildeqg() const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
     return 0.0;
  return Pqg() * (log((1.-z)/factory()->alphaParameter()) + log(min(1.,factory()->alphaParameter()/(1-z))));
}

double DipolePKOperator::KTildegq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return Pgq() * (log((1.-z)/factory()->alphaParameter())+ log(min(1.,factory()->alphaParameter()/(1-z))));
}




double DipolePKOperator::Pqq() const {
  assert(abs(parton->id()) < 7);
  double res = (3./2.+2.*log(1.-x)) * PDFx(parton);
  if ( z > x ) {
    res += 2.*(PDFxByz(parton) - z*PDFx(parton))/(z*(1.-z));
    res -= PDFxByz(parton) * (1.+z)/z;
  }
  return (CF*res);
}

double DipolePKOperator::Pqg() const {
  assert(parton->id() == ParticleID::g);
  if ( z < x )
    return 0.0;
  double res = 0.0;
  double factor = CF * ( 1. + sqr(1.-z) ) / sqr(z);
  // // int nlp = NLightProtonVec().size();
  // int nlp = lastBorn()->nLightProtonVec().size();
  // for ( int f = -nlp; f <= nlp; ++f ) {
  //   if ( f == 0 )
  //     continue;
  //   res += PDFxByz(getParticleData(f))*factor;
  // }
  for ( size_t f=0; f!=lastBorn()->nLightProtonVec().size(); ++f ) {
    res += PDFxByz(getParticleData(lastBorn()->nLightProtonVec()[f]))*factor;
    res += PDFxByz(getParticleData(-lastBorn()->nLightProtonVec()[f]))*factor;
  }
  return res;
}

double DipolePKOperator::Pgq() const {
  assert(abs(parton->id()) < 7);
  if ( z < x )
    return 0.0;
  return 0.5 * ( sqr(z) + sqr(1.-z) ) * PDFxByz(getParticleData(ParticleID::g)) / z;
}

double DipolePKOperator::Pgg() const {
  assert(parton->id() == ParticleID::g);
  double res = 
    // ( (11./6.) * CA - (1./3.)*NLightJetVec().size() + 2.*CA*log(1.-x) ) * PDFx(parton);
    ( (11./6.) * CA - (1./3.)*lastBorn()->nLightJetVec().size() + 2.*CA*log(1.-x) ) * PDFx(parton);
  if ( z > x ) {
    res += 2. * CA * ( PDFxByz(parton) - z*PDFx(parton) ) / (z*(1.-z));
    res += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

//////////////////////////////////////////////////////////////////////

// For every new momentum fraction at which we want to evaluate a pdf
// we introduce a new member to the pdf cache.

double DipolePKOperator::PDFx(tcPDPtr pd) const {
  map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator
    cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = make_pair(0.0,0.0);
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  // The convention of the pdf sets is always to return a*f(a). Upon usage, 
  // we have to remember that and rescale the result of the pdf accordingly
  // again, i.e. a*f(a)/a.
  if ( cached->second.first == 0.0 )
    cached->second.first = 
      pdf->xfx(particle,pd,scale,x)/x;
  return cached->second.first;
}

double DipolePKOperator::PDFxByz(tcPDPtr pd) const {
  map<pair<tcPDFPtr,tcPDPtr>,pair<double,double> >::iterator
    cached = pdfCache.find(make_pair(pdf,pd));
  if ( cached == pdfCache.end() ) {
    pdfCache[make_pair(pdf,pd)] = make_pair(0.0,0.0);
    cached = pdfCache.find(make_pair(pdf,pd));
  }
  // The convention of the pdf sets is always to return a*f(a). Upon usage, 
  // we have to remember that and rescale the result of the pdf accordingly
  // again, i.e. a*f(a)/a.
  if ( cached->second.second == 0.0 )
    cached->second.second = 
      pdf->xfx(particle,pd,scale,x/z)*z/x;
  return cached->second.second;
}

//////////////////////////////////////////////////////////////////////

void DipolePKOperator::persistentOutput(PersistentOStream & os) const {
  os << CA << CF << gammaQuark << gammaGluon << KQuark << KGluon
     << ounit(scale,GeV2) << pdf << particle << x << z << pdfCache << parton;
}

void DipolePKOperator::persistentInput(PersistentIStream & is, int) {
  is >> CA >> CF >> gammaQuark >> gammaGluon >> KQuark >> KGluon
     >> iunit(scale,GeV2) >> pdf >> particle >> x >> z >> pdfCache >> parton;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<DipolePKOperator,MatchboxInsertionOperator>
describeHerwigDipolePKOperator("Herwig::DipolePKOperator", "Herwig.so");

void DipolePKOperator::Init() {

  static ClassDocumentation<DipolePKOperator> documentation
    ("DipolePKOperator");

  DipoleRepository::registerInsertionPKOperator<0,DipolePKOperator>("LightPKOperator");

}

