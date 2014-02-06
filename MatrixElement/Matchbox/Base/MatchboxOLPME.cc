// -*- C++ -*-
//
// MatchboxOLPME.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxOLPME class.
//

#include "MatchboxOLPME.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "MatchboxMEBase.h"
#include "Herwig++/MatrixElement/Matchbox/MatchboxFactory.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxOLPME::MatchboxOLPME() 
  : theOrderInGs(0), theOrderInGem(0) {}

MatchboxOLPME::~MatchboxOLPME() {}

bool MatchboxOLPME::canHandle(const PDVector& p,
			      Ptr<MatchboxFactory>::tptr factory) const {

  if ( factory->processData()->diagramMap().find(p) !=
       factory->processData()->diagramMap().end() )
    return true;

  vector<Ptr<Tree2toNDiagram>::ptr> diags =
    factory->diagramGenerator()->generate(p,orderInGs(),orderInGem());

  if ( diags.empty() )
    return false;

  factory->processData()->diagramMap()[p] = diags;

  return true;

}

void MatchboxOLPME::setXComb(tStdXCombPtr xc) {
  theLastXComb = xc;
  lastMatchboxXComb(xc);
}

double MatchboxOLPME::me2() const {
  if ( !calculateTreeME2() )
    return crossingSign()*lastTreeME2();
  evalSubProcess();
  return crossingSign()*lastTreeME2();
}

double MatchboxOLPME::colourCorrelatedME2(pair<int,int> ij) const {
  double cfac = 1.;
  double Nc = generator()->standardModel()->Nc();
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  } else assert(false);
  if ( !calculateColourCorrelator(ij) )
    return crossingSign()*lastColourCorrelator(ij)/cfac;
  evalColourCorrelator(ij);
  return crossingSign()*lastColourCorrelator(ij)/cfac;
}

double MatchboxOLPME::spinColourCorrelatedME2(pair<int,int> ij,
					      const SpinCorrelationTensor& c) const {

  Lorentz5Momentum p = meMomenta()[ij.first];
  Lorentz5Momentum n = meMomenta()[ij.second];

  LorentzVector<Complex> polarization = plusPolarization(p,n,ij.first);

  Complex pFactor = (polarization*c.momentum())/sqrt(abs(c.scale()));

  double avg =
    colourCorrelatedME2(ij)*(-c.diagonal()+ (c.scale() > ZERO ? 1. : -1.)*norm(pFactor));

  Complex csCorr = 0.0;

  if ( calculateColourSpinCorrelator(ij) )
    evalSpinColourCorrelator(ij);

  csCorr = lastColourSpinCorrelator(ij);

  double corr = 
    2.*real(csCorr*sqr(pFactor));

  double Nc = generator()->standardModel()->Nc();
  double cfac = 1.;
  if ( mePartonData()[ij.first]->iColour() == PDT::Colour8 ) {
    cfac = Nc;
  } else if ( mePartonData()[ij.first]->iColour() == PDT::Colour3 ||
	      mePartonData()[ij.first]->iColour() == PDT::Colour3bar ) {
    cfac = (sqr(Nc)-1.)/(2.*Nc);
  } else assert(false);

  return 
    avg + crossingSign()*(c.scale() > ZERO ? 1. : -1.)*corr/cfac;

}

double MatchboxOLPME::oneLoopDoublePole() const {
  if ( !calculateOneLoopPoles() )
    return crossingSign()*lastOneLoopPoles().first;
  evalSubProcess();
  return crossingSign()*lastOneLoopPoles().first;
}

double MatchboxOLPME::oneLoopSinglePole() const {
  if ( !calculateOneLoopPoles() )
    return crossingSign()*lastOneLoopPoles().second;
  evalSubProcess();
  return crossingSign()*lastOneLoopPoles().second;
}

double MatchboxOLPME::oneLoopInterference() const {
  if ( !calculateOneLoopInterference() )
    return crossingSign()*lastOneLoopInterference();
  evalSubProcess();
  return crossingSign()*lastOneLoopInterference();
}

double MatchboxOLPME::largeNColourCorrelatedME2(pair<int,int>,
						Ptr<ColourBasis>::tptr) const {
  throw Exception() << "largeNColourCorrelatedME2 not supported by MatchboxOLPME"
		    << Exception::abortnow;
  return 0.;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

bool didstartOLP=false;

void MatchboxOLPME::doinitrun() {
  string contractFileName = name() + ".OLPContract.lh";
  int status = -1;
  if (didstartOLP==false){
	  startOLP(contractFileName,status);
	  didstartOLP=true;

  if ( status != 1 ) {
    throw Exception()
      << "Failed to restart one loop provider for amplitude '"
      << name() << "'\n" << Exception::abortnow;
  }
  }
  MatchboxAmplitude::doinitrun();
}

void MatchboxOLPME::persistentOutput(PersistentOStream & os) const {
  os << theOrderInGs << theOrderInGem;
}

void MatchboxOLPME::persistentInput(PersistentIStream & is, int) {
  is >> theOrderInGs >> theOrderInGem;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxOLPME,MatchboxAmplitude>
  describeHerwigMatchboxOLPME("Herwig::MatchboxOLPME", "HwMatchbox.so");

void MatchboxOLPME::Init() {

  static ClassDocumentation<MatchboxOLPME> documentation
    ("MatchboxOLPME implements OLP interfaces.");

}

