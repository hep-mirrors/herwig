// -*- C++ -*-
//
// QTildeFinder.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QTildeFinder class.
//

#include "QTildeFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/Event.h"
#include "Herwig/Shower/Base/ShowerParticle.h"
#include "ThePEG/Repository/UseRandom.h" 
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

DescribeClass<QTildeFinder,Herwig::PartnerFinder>
describeQTildeFinder ("Herwig::QTildeFinder","HwShower.so");

void QTildeFinder::persistentOutput(PersistentOStream & os) const {
  os << _finalFinalConditions << _initialFinalDecayConditions
     << _initialInitialConditions;
}

void QTildeFinder::persistentInput(PersistentIStream & is, int) {
  is >> _finalFinalConditions >> _initialFinalDecayConditions
     >>_initialInitialConditions;
}

void QTildeFinder::Init() {

  static ClassDocumentation<QTildeFinder> documentation
    ("This class is responsible for finding the partners for each interaction types ",
     "and within the evolution scale range specified by the ShowerVariables ",
     "then to determine the initial evolution scales for each pair of partners.");

  static Switch<QTildeFinder,unsigned int> interfaceFinalFinalConditions
    ("FinalFinalConditions",
     "The initial conditions for the shower of a final-final colour connection",
     &QTildeFinder::_finalFinalConditions, 0, false, false);
  static SwitchOption interfaceFinalFinalConditionsSymmetric
    (interfaceFinalFinalConditions,
     "Symmetric",
     "The symmetric choice",
     0);
  static SwitchOption interfaceFinalFinalConditionsColoured
    (interfaceFinalFinalConditions,
     "Coloured",
     "Maximal radiation from the coloured particle",
     1);
  static SwitchOption interfaceFinalFinalConditionsAntiColoured
    (interfaceFinalFinalConditions,
     "AntiColoured",
     "Maximal emission from the anticoloured particle",
     2);
  static SwitchOption interfaceFinalFinalConditionsRandom
    (interfaceFinalFinalConditions,
     "Random",
     "Randomly selected maximal emission from one of the particles",
     3);

  static Switch<QTildeFinder,unsigned int> interfaceInitialFinalDecayConditions
    ("InitialFinalDecayConditions",
     "The initial conditions for the shower of an initial-final"
     " decay colour connection.",
     &QTildeFinder::_initialFinalDecayConditions, 0, false, false);
  static SwitchOption interfaceInitialFinalDecayConditionsSymmetric
    (interfaceInitialFinalDecayConditions,
     "Symmetric",
     "The symmetric choice",
     0);
  static SwitchOption interfaceInitialFinalDecayConditionsMaximal
    (interfaceInitialFinalDecayConditions,
     "Maximal",
     "Maximal radiation from the decay product",
     1);
  static SwitchOption interfaceInitialFinalDecayConditionsSmooth
    (interfaceInitialFinalDecayConditions,
     "Smooth",
     "Smooth matching in the soft limit",
     2);

  static Switch<QTildeFinder,unsigned int> interfaceInitialInitialConditions
    ("InitialInitialConditions",
     "The initial conditions for the shower of an initial-initial"
     " colour connection.",
     &QTildeFinder::_initialInitialConditions, 0, false, false);
  static SwitchOption interfaceInitialInitialConditionsSymmetric
    (interfaceInitialInitialConditions,
     "Symmetric",
     "The symmetric choice",
     0);
  static SwitchOption interfaceInitialInitialConditionsMaximiseB
    (interfaceInitialInitialConditions,
     "MaximiseB",
     "Maximal radiation from parton b",
     1);
  static SwitchOption interfaceInitialInitialConditionsMaximiseC
    (interfaceInitialInitialConditions,
     "MaximiseC",
     "Maximal radiation from parton c",
     2);
}

pair<Energy,Energy> QTildeFinder::
calculateInitialFinalScales(const ShowerPPair &ppair, const bool isDecayCase) {
  return
    calculateInitialFinalScales(ppair.first->momentum(),ppair.second->momentum(),isDecayCase);
}

pair<Energy,Energy> QTildeFinder::
calculateInitialFinalScales(const Lorentz5Momentum& pb, const Lorentz5Momentum& pc,
			    const bool isDecayCase) {
  if(!isDecayCase) { 
    // In this case from JHEP 12(2003)045 we find the conditions
    // ktilde_b = (1+c) and ktilde_c = (1+2c)
    // We also find that c = m_c^2/Q^2. The process is a+b->c where
    // particle a is not colour connected (considered as a colour singlet).
    // Therefore we simply find that q_b = sqrt(Q^2+m_c^2) and 
    // q_c = sqrt(Q^2+2 m_c^2)
    // We also assume that the first particle in the pair is the initial
    // state particle and the second is the final state one c 
    Energy2  mc2 = sqr(pc.mass());
    Energy2  Q2  = -(pb-pc).m2();
    return pair<Energy,Energy>(sqrt(Q2+mc2), sqrt(Q2+2*mc2));
  }
  else {    
    // In this case from JHEP 12(2003)045 we find, for the decay
    // process b->c+a(neutral), the condition
    // (ktilde_b-1)*(ktilde_c-c)=(1/4)*sqr(1-a+c+lambda). 
    // We also assume that the first particle in the pair is the initial
    // state particle (b) and the second is the final state one (c).
    //  - We find maximal phase space coverage through emissions from 
    //    c if we set ktilde_c = 4.*(sqr(1.-sqrt(a))-c)
    //  - We find the most 'symmetric' way to populate the phase space
    //    occurs for (ktilde_b-1)=(ktilde_c-c)=(1/2)*(1-a+c+lambda) 
    //  - We find the most 'smooth' way to populate the phase space
    //    occurs for...
    Energy2 mb2(sqr(pb.mass()));
    double a=(pb-pc).m2()/mb2;
    double c=sqr(pc.mass())/mb2;
    double lambda   = 1. + a*a + c*c - 2.*a - 2.*c - 2.*a*c;
    lambda = sqrt(max(lambda,0.));
    double PROD     = 0.25*sqr(1. - a + c + lambda);
    double ktilde_b, ktilde_c,cosi(0.);
    switch(initialFinalDecayConditions()) {
    case 0: // the 'symmetric' choice
      ktilde_c = 0.5*(1-a+c+lambda) + c ;
      ktilde_b = 1.+PROD/(ktilde_c-c)   ;
      break;
    case 1:  // the 'maximal' choice
      ktilde_c = 4.0*(sqr(1.-sqrt(a))-c);
      ktilde_b = 1.+PROD/(ktilde_c-c)   ;
      break;
    case 2:  // the 'smooth' choice
      // c is a problem if very small here use 1GeV as minimum
      c = max(c,1.*GeV2/mb2);
      cosi = (sqr(1-sqrt(c))-a)/lambda;
      ktilde_b = 2.0/(1.0-cosi);
      ktilde_c = (1.0-a+c+lambda)*(1.0+c-a-lambda*cosi)/(2.0*(1.0+cosi));
      break;
    default:
      throw Exception() << "Invalid option for decay shower's phase space"
			<< " QTildeFinder::calculateInitialFinalScales"
			<< Exception::abortnow;
    }
    return pair<Energy,Energy>(sqrt(mb2*ktilde_b),sqrt(mb2*ktilde_c));
  }
}

pair<Energy,Energy> QTildeFinder::
calculateInitialInitialScales(const ShowerPPair &ppair) {
  return
    calculateInitialInitialScales(ppair.first->momentum(),
				  ppair.second->momentum());
}

pair<Energy,Energy> QTildeFinder::
calculateInitialInitialScales(const Lorentz5Momentum& p1, const Lorentz5Momentum& p2) {
  // This case is quite simple. From JHEP 12(2003)045 we find the condition
  // that ktilde_b = ktilde_c = 1. In this case we have the process
  // b+c->a so we need merely boost to the CM frame of the two incoming
  // particles and then qtilde is equal to the energy in that frame
  Energy Q = sqrt((p1+p2).m2());
  if(_initialInitialConditions==1) {
    return pair<Energy,Energy>(sqrt(2.0)*Q,sqrt(0.5)*Q);
  } else if(_initialInitialConditions==2) {
    return pair<Energy,Energy>(sqrt(0.5)*Q,sqrt(2.0)*Q);
  } else {
    return pair<Energy,Energy>(Q,Q);
  }
}

pair<Energy,Energy> QTildeFinder::
calculateFinalFinalScales(const ShowerPPair & pp) {
  bool colouredFirst =
    pp.first->colourLine()&&
    pp.first->colourLine()==pp.second->antiColourLine();
  return calculateFinalFinalScales(pp.first->momentum(),pp.second->momentum(),
				   colouredFirst);
}

pair<Energy,Energy> QTildeFinder::
calculateFinalFinalScales(Lorentz5Momentum p1, Lorentz5Momentum p2,
			  bool colouredFirst) {
  static const double eps=1e-7;
  // Using JHEP 12(2003)045 we find that we need ktilde = 1/2(1+b-c+lambda)
  // ktilde = qtilde^2/Q^2 therefore qtilde = sqrt(ktilde*Q^2)
  // find momenta in rest frame of system
  // calculate quantities for the scales
  Energy2 Q2 = (p1+p2).m2();
  double b = p1.mass2()/Q2;
  double c = p2.mass2()/Q2;
  if(b<0.) {
    if(b<-eps) {
      throw Exception() << "Negative Mass squared b = " << b
			<< "in QTildeFinder::calculateFinalFinalScales()"
			<< Exception::eventerror;
    }
    b = 0.;
  }
  if(c<0.) {
    if(c<-eps) {
      throw Exception() << "Negative Mass squared c = " << c
			<< "in QTildeFinder::calculateFinalFinalScales()"
			<< Exception::eventerror;
    }
    c = 0.;
  }
  // KMH & PR - 16 May 2008 - swapped lambda calculation from 
  // double lam=2.*p1.vect().mag()/Q; to sqrt(kallen(1,b,c)), 
  // which should be identical for p1 & p2 onshell in their COM
  // but in the inverse construction for the Nason method, this
  // was not the case, leading to misuse. 
  double lam=sqrt((1.+sqrt(b)+sqrt(c))*(1.-sqrt(b)-sqrt(c))
                 *(sqrt(b)-1.-sqrt(c))*(sqrt(c)-1.-sqrt(b)));
  // symmetric case
  unsigned int iopt=finalFinalConditions();
  Energy firstQ,secondQ;
  if(iopt==0) {
    firstQ  = sqrt(0.5*Q2*(1.+b-c+lam));
    secondQ = sqrt(0.5*Q2*(1.-b+c+lam));
  }
  // assymetric choice
  else {
    double kappab,kappac;
    // calculate kappa with coloured line getting maximum
    if((iopt==1&&colouredFirst)|| // first particle coloured+maximal for coloured
       (iopt==2&&!colouredFirst)|| // first particle anticoloured+maximal for acoloured
       (iopt==3&&UseRandom::rndbool(0.5))) { // random choice
      kappab=4.*(1.-2.*sqrt(c)-b+c);
      kappac=c+0.25*sqr(1.-b-c+lam)/(kappab-b);
    }
    else {
      kappac=4.*(1.-2.*sqrt(b)-c+b);
      kappab=b+0.25*sqr(1.-b-c+lam)/(kappac-c);
    }
    // calculate the scales
    firstQ  = sqrt(Q2*kappab);
    secondQ = sqrt(Q2*kappac);
  }
  return pair<Energy,Energy>(firstQ, secondQ);
}
