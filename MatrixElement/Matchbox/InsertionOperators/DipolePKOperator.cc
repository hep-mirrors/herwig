// -*- C++ -*-
//
// DipolePKOperator.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;
using Constants::pi;

DipolePKOperator::DipolePKOperator() 
  : MatchboxInsertionOperator(),
    CA(-1.0), CF(-1.0), 
    gammaQuark(-1.0), gammaGluon(-1.0),
    KQuark(-1.0), KGluon(-1.0) {}

DipolePKOperator::~DipolePKOperator() {}

IBPtr DipolePKOperator::clone() const {
  return new_ptr(*this);
}

IBPtr DipolePKOperator::fullclone() const {
  return new_ptr(*this);
}

void DipolePKOperator::setXComb(tStdXCombPtr xc) {
  MatchboxInsertionOperator::setXComb(xc);
  if ( CA < 0. ) {
    CA = SM().Nc();
    CF = (SM().Nc()*SM().Nc()-1.0)/(2.*SM().Nc());
    gammaQuark = (3./2.)*CF;
    gammaGluon = (11./6.)*CA - (1./3.)*lastBorn()->nLight();
    KQuark = (7./2.-sqr(pi)/6.)*CF;
    KGluon = (67./18.-sqr(pi)/6.)*CA-(5./9.)*lastBorn()->nLight();
  }
}

double DipolePKOperator::gammaSoft() const {
  double res = (1.+log(1.-x))*PDFx(parton);
  if ( z > x )
    res +=
      (PDFxByz(parton) - z*PDFx(parton)) / (z*(1.-z));
  return res;
}

double DipolePKOperator::softLogByz(tcPDPtr p) const {
  double res = ( sqr(log(1.-x))/2. - sqr(pi)/6. ) * PDFx(p);
  if ( z > x ) {
    res += (PDFxByz(p) - z*PDFx(p))*log(1.-z)/(z*(1.-z));
    res -= PDFxByz(p)*log(z)/(z*(1.-z));
  }
  return res;
}

double DipolePKOperator::softLog(tcPDPtr p) const {
  double res = sqr(log(1.-x)) * PDFx(p) / 2.;
  if ( z > x ) {
    res += (PDFxByz(p) - z*PDFx(p))*log(1.-z)/(z*(1.-z));
  }
  return res;
}

double DipolePKOperator::KBarqq() const {
  assert(abs(parton->id()) < 6);
  double res = 
    2.*softLogByz(parton) +
    (sqr(pi) - 5.)*PDFx(parton);
  if ( z > x ) {
    res += PDFxByz(parton)*( (1.-z) - (1.+z)*log((1.-z)/z) ) / z;
  }
  return (res * CF);
}

double DipolePKOperator::KTildeqq() const {
  assert(abs(parton->id()) < 6);
  double res =
    2.*CF*softLog(parton) - CF*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res -= ( CF * (1.+z) * log(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

double DipolePKOperator::Pqq() const {
  assert(abs(parton->id()) < 6);
  double res = (3./2.+2.*log(1.-x)) * PDFx(parton);
  if ( z > x ) {
    res += 2.*(PDFxByz(parton) - z*PDFx(parton))/(z*(1.-z));
    res -= PDFxByz(parton) * (1.+z)/z;
  }
  return (CF*res);
}

double DipolePKOperator::KBarqg() const {
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

double DipolePKOperator::KTildeqg() const {
  if ( z < x )
    return 0.0;
  return Pqg() * log(1.-z);
}

double DipolePKOperator::Pqg() const {
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

double DipolePKOperator::KBargq() const {
  assert(abs(parton->id()) < 6);
  if ( z < x )
    return 0.0;
  return
    PDFxByz(getParticleData(ParticleID::g)) *
    ( 0.5*(sqr(z)+sqr(1.-z))*log((1.-z)/z) + z*(1.-z) ) / z;
}

double DipolePKOperator::KTildegq() const {
  assert(abs(parton->id()) < 6);
  if ( z < x )
    return 0.0;
  return Pgq() * log(1.-z);
}

double DipolePKOperator::Pgq() const {
  assert(abs(parton->id()) < 6);
  if ( z < x )
    return 0.0;
  return 0.5 * ( sqr(z) + sqr(1.-z) ) * PDFxByz(getParticleData(ParticleID::g)) / z;
}

double DipolePKOperator::KBargg() const {
  assert(parton->id() == ParticleID::g);
  double res = 
    2.* CA* softLogByz(parton) +
    ( CA*( sqr(pi) - 50./9. ) + (8./9.)*lastBorn()->nLight() ) * PDFx(parton);
  if ( z > x ) {
    res += 2.*CA*((1.-z)/z-1.+z*(1.-z))*log((1.-z)/z)*PDFxByz(parton)/z;
  }
  return res;
}

double DipolePKOperator::KTildegg() const {
  assert(parton->id() == ParticleID::g);
  double res =
    2.*CA*softLog(parton) - CA*(sqr(pi)/3.)*PDFx(parton);
  if ( z > x ) {
    res += ( 2.*CA * ( (1.-z)/z -1. + z*(1.-z) ) * log(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

double DipolePKOperator::Pgg() const {
  assert(parton->id() == ParticleID::g);
  double res = 
    ( (11./6.) * CA - (1./3.) * lastBorn()->nLight() + 2.*CA*log(1.-x) ) * PDFx(parton);
  if ( z > x ) {
    res += 2. * CA * ( PDFxByz(parton) - z*PDFx(parton) ) / (z*(1.-z));
    res += 2.* CA *( (1.-z)/z - 1. + z*(1.-z) ) * PDFxByz(parton) / z;
  }
  return res;
}

double DipolePKOperator::PDFx(tcPDPtr pd) const {
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

double DipolePKOperator::PDFxByz(tcPDPtr pd) const {
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

bool DipolePKOperator::apply(tcPDPtr pd) const {
  return
    pd->mass() == ZERO &&
    (abs(pd->id()) < 6 || pd->id() == ParticleID::g);
}

bool DipolePKOperator::apply(const cPDVector& pd) const {
  if ( !apply(pd[0]) && !apply(pd[1]) )
    return false;
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

double DipolePKOperator::sumParton(int id) const {

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

  if ( mePartonData()[id]->id() == ParticleID::g )
    res += (KBargg() + KBarqg())*lastBorn()->me2();

  if ( abs(mePartonData()[id]->id()) < 6 )
    res += (KBarqq() + KBargq())*lastBorn()->me2();

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

    if ( abs(mePartonData()[id]->id()) < 6 &&
	 thePqq == 0.0 ) {
      thePgq = Pgq();
      thePqq = Pqq();
    }

    double theLog = log(scale/(2.*((*Pi)*meMomenta()[id])));

    res +=
      ( (**i).id() == ParticleID::g ? gammaGluon : gammaQuark ) *
      theGammaSoft * fiCorrelated;

    if ( mePartonData()[id]->id() == ParticleID::g ) {
      res +=
	( thePgg + thePqg ) * theLog * ifCorrelated;
    }

    if ( abs(mePartonData()[id]->id()) < 6 ) {
      res +=
	( thePqq + thePgq ) * theLog * ifCorrelated;
      if ( disFinite == 0.0 && z > x ) {
	disFinite = CF*PDFxByz(parton)*(1.+3.*z/2.)/z;
      }
      if ( z > x )
	res -= disFinite*ifCorrelated;
    }

    if ( abs(mePartonData()[idi]->id()) < 6 ) {
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

  }

  if ( mePartonData()[ id == 0 ? 1 : 0 ]->coloured() &&
       !lastBorn()->noDipole(id == 0 ? 0 : 1,
			     id == 0 ? 1 : 0) ) {

    if ( mePartonData()[id]->id() == ParticleID::g &&
	 thePgg == 0.0 ) {
      thePqg = Pqg();
      thePgg = Pgg();
    }

    if ( abs(mePartonData()[id]->id()) < 6 &&
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

    if ( abs(mePartonData()[id]->id()) < 6 ) {
      res +=
	( thePqq + thePgq ) * theLog * iiCorrelated;
      res -=
	( KTildeqq() + KTildegq() ) * iiCorrelated;
    }

  }

  return res * mapz;

}

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

  DipoleRepository::registerInsertionOperator<0,DipolePKOperator>("LightPKOperator");

}

