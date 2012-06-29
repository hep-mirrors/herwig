// -*- C++ -*-
//
// MatchboxMEllbar2qqbarg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMEllbar2qqbarg class.
//

#include "MatchboxMEllbar2qqbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

MatchboxMEllbar2qqbarg::MatchboxMEllbar2qqbarg() 
  : MatchboxMEBase(), theUserScale(0.0*GeV) {}

MatchboxMEllbar2qqbarg::~MatchboxMEllbar2qqbarg() {}

IBPtr MatchboxMEllbar2qqbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMEllbar2qqbarg::fullclone() const {
  return new_ptr(*this);
}


bool MatchboxMEllbar2qqbarg::generateKinematics(const double * r) {

  if ( phasespace() )
    return MatchboxMEBase::generateKinematics(r);

  using namespace RandomHelpers;

  double x0 = .01;
  double xc = 1e-4;

  pair<double,double> xw1 = generate((piecewise(),
				      inverse(1.,0.,1.-x0),
				      match(flat(1.-x0,1.-xc))),
				     r[0]);
  double x1 = xw1.first;
  pair<double,double> xw2 = generate((piecewise(),
				      inverse(1.,0.,1.-x0),
				      match(flat(1.-x0,1.-xc))),
				     r[1]);
  double x2 = xw2.first;

  double mapX = xw1.second*xw2.second;
  
  double mu2 = sqr( getParticleData( mePartonData()[2]->id() )->mass() ) / lastSHat();

  if ( x1 < 2.*sqrt(mu2) || x1 > 1. ||
    x2 < ( (2.-x1)*(1.-x1+2.*mu2) - (1.-x1)*sqrt(x1*x1-4.*mu2) ) / ( 2.*(1.-x1+mu2) ) ||
    x2 > ( (2.-x1)*(1.-x1+2.*mu2) + (1.-x1)*sqrt(x1*x1-4.*mu2) ) / ( 2.*(1.-x1+mu2) ) ) {
    jacobian(0.0);
    return false;
  }

  double x3 = 2.-x1-x2;
  
  double th1 = acos(2.*r[2]-1.);
  double phi1 = 2.*Constants::pi*r[3];
  double phi2 = 2.*Constants::pi*r[4];
  
  double ct12 = (sqr(x3)-sqr(x1)-sqr(x2)+8.*mu2) /
    ( 2. * sqrt(sqr(x1)-4.*mu2) * sqrt(sqr(x2)-4.*mu2) );
  double st12 = sqrt(1.-sqr(ct12));

  Axis n1 (0.,0.,1.);
  Axis n2 (st12,0.,ct12);
  Axis x3n3 = -sqrt(sqr(x1)-4.*mu2)*n1 -sqrt(sqr(x2)-4.*mu2)*n2;

  n2.rotate(phi1,n1);
  x3n3.rotate(phi1,n1);

  Axis xPrime (1.,0.,0.);
  xPrime.rotate(phi1,n1);

  n1.rotate(th1,xPrime);
  n2.rotate(th1,xPrime);
  x3n3.rotate(th1,xPrime);

  Axis zPrime = n1;

  n1.rotate(phi2,zPrime);
  n2.rotate(phi2,zPrime);
  x3n3.rotate(phi2,zPrime);

  meMomenta()[2] = Lorentz5Momentum(sqrt(mu2*lastSHat()),(sqrt(lastSHat())*sqrt(sqr(x1)-4.*mu2)/2.)*n1);
  meMomenta()[3] = Lorentz5Momentum(sqrt(mu2*lastSHat()),(sqrt(lastSHat())*sqrt(sqr(x2)-4.*mu2)/2.)*n2);
  meMomenta()[4] = Lorentz5Momentum(ZERO,(sqrt(lastSHat())/2.)*x3n3);

  jacobian(mapX/(128.*Constants::pi*Constants::pi*Constants::pi));

  setScale();
  logGenerateKinematics(r);
  return true;

}

double MatchboxMEllbar2qqbarg::me2() const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::me2();

  double res;
  if ( !calculateME2(res) )
    return res;

  Lorentz5Momentum pq = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[2] : 
    meMomenta()[3];

  Lorentz5Momentum pqbar = 
    mePartonData()[2]->id() > 0 ? 
    meMomenta()[3] : 
    meMomenta()[2];

  Lorentz5Momentum pl = 
    mePartonData()[0]->id() > 0 ? 
    - meMomenta()[1] : 
    - meMomenta()[0];

  Lorentz5Momentum plbar = 
    mePartonData()[0]->id() > 0 ? 
    - meMomenta()[0] : 
    - meMomenta()[1];

  Lorentz5Momentum pg = meMomenta()[4];

  prepare(pl,plbar,pq,pqbar,pg,lastSHat(),mePartonData()[0],mePartonData()[2]);

  lastME2(MatchboxMEllbarqqbarg::evaluateME2()*me2Norm());
  cacheME2(lastME2());

  logME2();

  return lastME2();

}

Energy2 MatchboxMEllbar2qqbarg::factorizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

Energy2 MatchboxMEllbar2qqbarg::renormalizationScale() const {

  if ( theUserScale != ZERO )
    return sqr(theUserScale);

  return lastSHat();

}

double MatchboxMEllbar2qqbarg::colourCorrelatedME2(pair<int,int> ij) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::colourCorrelatedME2(ij);

  generator()->logWarning(Exception() 
			  << "A non-exisiting colour correlation was requested "
			  << "from the matrix element '" << name() << "'."
			  << Exception::warning);
  lastME2(0.0);
  return lastME2();

}

double MatchboxMEllbar2qqbarg::spinColourCorrelatedME2(pair<int,int> ij,
						       const SpinCorrelationTensor& c) const {

  if ( matchboxAmplitude() )
    return MatchboxMEBase::spinColourCorrelatedME2(ij,c);

  generator()->logWarning(Exception() 
			  << "A non-exisiting colour correlation was requested "
			  << "from the matrix element '" << name() << "'."
			  << Exception::warning);

  lastME2(0.0);
  return lastME2();

}

void MatchboxMEllbar2qqbarg::getDiagrams() const {

  tcPDPtr g = getParticleData(ParticleID::g);
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);

  for ( PDVector::const_iterator l = theLeptonFlavours.begin();
	l != theLeptonFlavours.end(); ++l )
    for ( PDVector::const_iterator q = theQuarkFlavours.begin();
	  q != theQuarkFlavours.end(); ++q ) {
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, gamma, 3, *q, 4, *q, 3, (**q).CC(), 4, g, -1)));
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, gamma, 3, *q, 3, (**q).CC(), 5, (**q).CC(), 5, g, -2)));
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, Z0, 3, *q, 4, *q, 3, (**q).CC(), 4, g, -3)));
      add(new_ptr((Tree2toNDiagram(2), *l, (**l).CC(), 1, Z0, 3, *q, 3, (**q).CC(), 5, (**q).CC(), 5, g, -4)));
    }

}

Selector<MEBase::DiagramIndex> 
MatchboxMEllbar2qqbarg::diagrams(const DiagramVector &) const {
  Selector<MEBase::DiagramIndex> sel;
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  double wGamma = 
    sqr( sqr(generator()->maximumCMEnergy())/lastSHat() );
  double wZ = 
    sqr( sqr(generator()->maximumCMEnergy()) ) /
    ( sqr(lastSHat()-sqr(Z0->mass())) + sqr(Z0->mass())*sqr(Z0->width()) );
  if ( meMomenta()[2]*meMomenta()[4] < meMomenta()[3]*meMomenta()[4] ) {
    sel.insert(wGamma,0);
    sel.insert(wZ,2);
  } else {
    sel.insert(wGamma,1);
    sel.insert(wZ,3);
  }
  return sel;
}

Selector<const ColourLines *>
MatchboxMEllbar2qqbarg::colourGeometries(tcDiagPtr d) const {
  static const ColourLines cq("5 -7, 7 4 -6"); 
  static const ColourLines cqbar("4 -5 -7, 7 -6"); 
  static const ColourLines cbarq("-5 7, -7 -4 6"); 
  static const ColourLines cbarqbar("-4 5 7, -7 6"); 
  Selector<const ColourLines *> sel;
  if ( mePartonData()[2]->id() > 0 ) {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cq);
    else
      sel.insert(1.0, &cqbar);
  } else {
    if ( abs(d->id()) % 2 != 0 )
      sel.insert(1.0, &cbarq);
    else
      sel.insert(1.0, &cbarqbar);
  }
  return sel;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxMEllbar2qqbarg::doinit() {
  MatchboxMEBase::doinit();
  MatchboxMEllbarqqbarg::doinit(*this);
  // comment out if working with massive quarks (flagMS1, also in FFM??xDipole.cc)
//   for ( PDVector::const_iterator q = theQuarkFlavours.begin();
// 	q != theQuarkFlavours.end(); ++q )
//     if ( (**q).mass() != ZERO )
//       Throw<InitException>() << "The matrix element '"
// 			     << name() << "' is only capable of "
// 			     << "producing massless quarks.";
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxMEllbar2qqbarg::persistentOutput(PersistentOStream & os) const {
  MatchboxMEllbarqqbarg::persistentOutput(os);
  os << theLeptonFlavours << theQuarkFlavours << ounit(theUserScale,GeV);
}

void MatchboxMEllbar2qqbarg::persistentInput(PersistentIStream & is, int) {
  MatchboxMEllbarqqbarg::persistentInput(is);
  is >> theLeptonFlavours >> theQuarkFlavours >> iunit(theUserScale,GeV);
}

ClassDescription<MatchboxMEllbar2qqbarg> MatchboxMEllbar2qqbarg::initMatchboxMEllbar2qqbarg;
// Definition of the static class description member.

void MatchboxMEllbar2qqbarg::Init() {

  static ClassDocumentation<MatchboxMEllbar2qqbarg> documentation
    ("MatchboxMEllbar2qqbarg");

  static RefVector<MatchboxMEllbar2qqbarg,ParticleData> interfaceLeptonFlavours
    ("LeptonFlavours",
     "The lepton flavours for this matrix element.",
     &MatchboxMEllbar2qqbarg::theLeptonFlavours, -1, false, false, true, true, false);

  static RefVector<MatchboxMEllbar2qqbarg,ParticleData> interfaceQuarkFlavours
    ("QuarkFlavours",
     "The quark flavours for this matrix element.",
     &MatchboxMEllbar2qqbarg::theQuarkFlavours, -1, false, false, true, true, false);


  static Parameter<MatchboxMEllbar2qqbarg,Energy> interfaceUserScale
    ("UserScale",
     "A user defined renormalization scale.",
     &MatchboxMEllbar2qqbarg::theUserScale, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

}

