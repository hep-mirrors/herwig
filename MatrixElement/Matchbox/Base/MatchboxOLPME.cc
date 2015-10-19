// -*- C++ -*-
//
// MatchboxOLPME.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxOLPME::MatchboxOLPME() 
  : theOrderInGs(0), theOrderInGem(0), theSetMuToMuR(false), 
    theUseRunningAlphaS(false), theUseRunningAlphaEW(false) {}

MatchboxOLPME::~MatchboxOLPME() {}

bool MatchboxOLPME::canHandle(const PDVector& p,
			      Ptr<MatchboxFactory>::tptr factory,
			      bool) const {

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
    return lastTreeME2();
  evalSubProcess();
  return lastTreeME2();
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
    return lastColourCorrelator(ij)/cfac;
  evalColourCorrelator(ij);
  return lastColourCorrelator(ij)/cfac;
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
    avg + (c.scale() > ZERO ? 1. : -1.)*corr/cfac;

}

double MatchboxOLPME::spinCorrelatedME2(pair<int,int> ij,
					      const SpinCorrelationTensor& c) const {

  Lorentz5Momentum p = meMomenta()[ij.first];
  Lorentz5Momentum n = meMomenta()[ij.second];

  LorentzVector<Complex> polarization = plusPolarization(p,n,ij.first);

  Complex pFactor = (polarization*c.momentum())/sqrt(abs(c.scale()));

  double avg =
    me2()*(-c.diagonal()+ (c.scale() > ZERO ? 1. : -1.)*norm(pFactor));

  Complex csCorr = 0.0;

  if ( calculateSpinCorrelator(ij) )
    evalSpinCorrelator(ij);

  csCorr = lastSpinCorrelator(ij);

  double corr = 
    2.*real(csCorr*sqr(pFactor));

  return 
    avg + (c.scale() > ZERO ? 1. : -1.)*corr;

}


void MatchboxOLPME::evalSpinCorrelator(pair<int,int>) const {
  throw Exception()
    << "MatchboxOLPME::spinCorrelatedME2() is not implemented.\n"
    << "Please check your setup." << Exception::runerror;
}

double MatchboxOLPME::oneLoopDoublePole() const {
  if ( !calculateOneLoopPoles() )
    return lastOneLoopPoles().first;
  evalSubProcess();
  return lastOneLoopPoles().first;
}

double MatchboxOLPME::oneLoopSinglePole() const {
  if ( !calculateOneLoopPoles() )
    return lastOneLoopPoles().second;
  evalSubProcess();
  return lastOneLoopPoles().second;
}

double MatchboxOLPME::oneLoopInterference() const {
  if ( !calculateOneLoopInterference() )
    return lastOneLoopInterference();
  evalSubProcess();
  return lastOneLoopInterference();
}

double MatchboxOLPME::largeNColourCorrelatedME2(pair<int,int> ij,
						Ptr<ColourBasis>::tptr basis) const {
  if ( trivialColourLegs() )
    return MatchboxAmplitude::largeNColourCorrelatedME2(ij,basis);
  throw Exception() << "MatchboxOLPME::largeNColourCorrelatedME2(): not supported"
		    << Exception::runerror;
  return 0.;
}

double MatchboxOLPME::largeNME2(Ptr<ColourBasis>::tptr basis) const {
  if ( trivialColourLegs() )
    return MatchboxAmplitude::largeNME2(basis);
  throw Exception() << "MatchboxOLPME::largeNME2(): not supported"
		    << Exception::runerror;
  return 0.;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxOLPME::doinit() {
  if ( theUseRunningAlphaS && !theSetMuToMuR ) {
    throw Exception() << "MatchboxOLPME::doinit(): "
    << "Amplitude '" << name() << "' "
    << "uses a running alpha_s but fixed renormalization scale!\n"
    << Exception::runerror;
  }
  if ( !theUseRunningAlphaS && theSetMuToMuR ) {
    throw Exception() << "MatchboxOLPME::doinit(): "
    << "Amplitude '" << name() << "' "
    << "uses a fixed alpha_s but running renormalization scale!\n"
    << Exception::runerror;
  }
  if ( !didStartOLP() ) {
    string contractFileName = 
      optionalContractFile().empty() ? 
      factory()->buildStorage() + name() + ".OLPContract.lh" :
      optionalContractFile();
    int status = -1;
    startOLP(contractFileName,status);
    didStartOLP()=true;
    if ( status != 1 ) {
      throw Exception() << "MatchboxOLPME::doinit(): "
	<< "Failed to restart one loop provider for amplitude '"
	<< name() << "'\n" << Exception::runerror;
    }
  }
  MatchboxAmplitude::doinit();
}

void MatchboxOLPME::doinitrun() {
  if ( theUseRunningAlphaS && !theSetMuToMuR ) {
    throw Exception() << "MatchboxOLPME::doinitrun(): "
    << "Amplitude '" << name() << "' "
    << "uses a running alpha_s but fixed renormalization scale!\n"
    << Exception::runerror;
  }
  if ( !theUseRunningAlphaS && theSetMuToMuR ) {
    throw Exception() << "MatchboxOLPME::doinitrun(): "
    << "Amplitude '" << name() << "' "
    << "uses a fixed alpha_s but running renormalization scale!\n"
    << Exception::runerror;
  }
  if ( !didStartOLP() ) {
    string contractFileName = 
      optionalContractFile().empty() ? 
      factory()->buildStorage() + name() + ".OLPContract.lh" :
      optionalContractFile();
    int status = -1;
    startOLP(contractFileName,status);
    didStartOLP()=true;
    if ( status != 1 ) {
      throw Exception() << "MatchboxOLPME::doinitrun(): "
	<< "Failed to restart one loop provider for amplitude '"
	<< name() << "'\n" << Exception::runerror;
    }
  }
  MatchboxAmplitude::doinitrun();
}


Energy2 MatchboxOLPME::mu2() const { 
  if (theSetMuToMuR) {
    return lastMatchboxXComb()->lastRenormalizationScale();
  }
  return lastSHat(); 
}

bool MatchboxOLPME::hasRunningAlphaS() const { 
  if (theUseRunningAlphaS) {
    return true;
  }
  return false; 
}

bool MatchboxOLPME::hasRunningAlphaEW() const { 
  if (theUseRunningAlphaEW) {
    return true;
  }
  return false; 
}

void MatchboxOLPME::persistentOutput(PersistentOStream & os) const {
  os << theOrderInGs << theOrderInGem << theSetMuToMuR << theUseRunningAlphaS << theUseRunningAlphaEW;
}

void MatchboxOLPME::persistentInput(PersistentIStream & is, int) {
  is >> theOrderInGs >> theOrderInGem >> theSetMuToMuR >> theUseRunningAlphaS >> theUseRunningAlphaEW;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxOLPME,MatchboxAmplitude>
  describeHerwigMatchboxOLPME("Herwig::MatchboxOLPME", "Herwig.so");

void MatchboxOLPME::Init() {

  static ClassDocumentation<MatchboxOLPME> documentation
    ("MatchboxOLPME implements OLP interfaces.");

  static Switch<MatchboxOLPME,bool> interfaceSetMuToMuR
         ("SetMuToMuR",
          "Switch On to set the value of the dimensional regularization parameter mu2 for this OLP"
          "to the value of the renormalization scale muR2. Default is Off. The restoration for the "
          "full renormalization scale dependence in the DipoleIOperator isn't needed in this case.",
          &MatchboxOLPME::theSetMuToMuR, false, false, false);
  static SwitchOption interfaceSetMuToMuROn
         (interfaceSetMuToMuR,
          "On",
          "On",
          true);
  static SwitchOption interfaceSetMuToMuROff
         (interfaceSetMuToMuR,
          "Off",
          "Off",
          false);

  static Switch<MatchboxOLPME,bool> interfaceUseRunningAlphaS
         ("UseRunningAlphaS",
          "Switch On to set the value of alpha_s for this OLP to the value of the running alpha_s "
          "instead of to the value of the reference alpha_s. Default is Off. This also sets the value "
          "for hasRunningAlphaS() to true.",
          &MatchboxOLPME::theUseRunningAlphaS, false, false, false);
  static SwitchOption interfaceUseRunningAlphaSOn
         (interfaceUseRunningAlphaS,
          "On",
          "On",
          true);
  static SwitchOption interfaceUseRunningAlphaSOff
         (interfaceUseRunningAlphaS,
          "Off",
          "Off",
          false);

  static Switch<MatchboxOLPME,bool> interfaceUseRunningAlphaEW
         ("UseRunningAlphaEW",
          "Switch On to set the value of alpha_ew for this OLP to the value of the running alpha_ew "
          "instead of to the value of the reference alpha_ew. Default is Off. This also sets the value "
          "for hasRunningAlphaEW() to true.",
          &MatchboxOLPME::theUseRunningAlphaEW, false, false, false);
  static SwitchOption interfaceUseRunningAlphaEWOn
         (interfaceUseRunningAlphaEW,
          "On",
          "On",
          true);
  static SwitchOption interfaceUseRunningAlphaEWOff
         (interfaceUseRunningAlphaEW,
          "Off",
          "Off",
          false);

}

