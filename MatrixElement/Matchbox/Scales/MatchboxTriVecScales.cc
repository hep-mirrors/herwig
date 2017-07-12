// -*- C++ -*-
//
// MatchboxTriVecScales.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxTriVecScales class.
//

#include "MatchboxTriVecScales.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxTriVecScales::MatchboxTriVecScales()
  : theTriVecScaleChoice(1) {}

MatchboxTriVecScales::~MatchboxTriVecScales() {}

IBPtr MatchboxTriVecScales::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxTriVecScales::fullclone() const {
  return new_ptr(*this);
}

Energy2 MatchboxTriVecScales::renormalizationScale() const {

  // calculate scales
  tcPDVector pd (mePartonData().begin() + 2, mePartonData().end());
  vector<LorentzMomentum> p (meMomenta().begin() + 2, meMomenta().end());
  tcPDPtr t1 = mePartonData()[0];
  tcPDPtr t2 = mePartonData()[1];
  tcCutsPtr cuts = lastCutsPtr();

  theJetFinder->cluster(pd, p, cuts, t1, t2);

  Energy sumpartpt = ZERO;
  Energy vecet[3] = {ZERO,ZERO,ZERO};
  Energy sumvecet = ZERO;
  int foundlept[3] = {0,0,0}; // First entry for e-like, second entry for mu-like, third entry for tau-like
  int ivec[3] = {0,0,0};
  tcPDVector::const_iterator itpd = pd.begin();
  int ip = 2;

  double avgy12 = 0;
  if ( theTriVecScaleChoice==2 ) {
    LorentzMomentum p1 = LorentzMomentum();
    LorentzMomentum p2 = LorentzMomentum();
    int njets = 0;
    for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
         itp != p.end(); ++itp, ++itpd, ++ip ) 
      if ( (**itpd).coloured() ) {
	++njets;
        if ( itp->perp() > p1.perp() ) {
	  p1 = *itp;
	  p2 = p1;
	} else if ( itp->perp() > p2.perp() ) 
	  p2 = *itp;
      }
    if ( njets < 2 )
      throw Exception() << "MatchboxTriVecScales: Not enough jets in event for HtPrimeModScale!"
			<< Exception::runerror;
    avgy12 = (p1.rapidity()+p2.rapidity())/2.;
  }

  for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
       itp != p.end(); ++itp, ++itpd, ++ip ) {
    if ( (**itpd).coloured() ) {
      if ( theTriVecScaleChoice==2 ) 
        sumpartpt += (*itp).perp()*exp(abs((*itp).rapidity()-avgy12));
      else
        sumpartpt += (*itp).perp();
    }
    if ( !(**itpd).coloured() ) {
      if ( abs((**itpd).id())==ParticleID::nu_e || abs((**itpd).id())==ParticleID::eminus ) {
        if ( ++foundlept[0] == 1 )
          ivec[0] = ip;
	if ( foundlept[0] == 2 )
          vecet[0] = sqrt( (p[ivec[0]]+p[ip]).m2() + (p[ivec[0]]+p[ip]).perp2() );
      }
      if ( abs((**itpd).id())==ParticleID::nu_mu || abs((**itpd).id())==ParticleID::muminus ) {
        if ( ++foundlept[1] == 1 )
          ivec[1] = ip;
	if ( foundlept[1] == 2 )
          vecet[1] = sqrt( (p[ivec[1]]+p[ip]).m2() + (p[ivec[1]]+p[ip]).perp2() );
      }
      if ( abs((**itpd).id())==ParticleID::nu_tau || abs((**itpd).id())==ParticleID::tauminus ) {
        if ( ++foundlept[2] == 1 )
          ivec[2] = ip;
	if ( foundlept[2] == 2 )
          vecet[2] = sqrt( (p[ivec[2]]+p[ip]).m2() + (p[ivec[2]]+p[ip]).perp2() );
      }
    }
  }

  // Check for consistency in number of lepton pairs and members therein
  for (int i = 0; i<3; ++i ) {
    if ( foundlept[i] != 2 ) 
      throw Exception() << "MatchboxTriVecScales: Inconsistency in number of lepton pairs and members therein!"
			<< Exception::runerror;
  }

  sumvecet = vecet[0]+vecet[1]+vecet[2];

  if ( theTriVecScaleChoice==1 || theTriVecScaleChoice==2 ) 
    return sqr(0.5*(sumpartpt+sumvecet));
  else if ( theTriVecScaleChoice==3 ) 
    return sqr(0.5*sumvecet);
  else
    throw Exception() << "MatchboxTriVecScales: Scale choice out of range: " 
	              << theTriVecScaleChoice
                      << Exception::runerror;
}

Energy2 MatchboxTriVecScales::factorizationScale() const {
  return renormalizationScale();
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxTriVecScales::persistentOutput(PersistentOStream & os) const {
  os << theJetFinder << theTriVecScaleChoice;
}

void MatchboxTriVecScales::persistentInput(PersistentIStream & is, int) {
  is >> theJetFinder >> theTriVecScaleChoice;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxTriVecScales,MatchboxScaleChoice>
  describeHerwigMatchboxTriVecScales("Herwig::MatchboxTriVecScales", "HwMatchboxScales.so");

void MatchboxTriVecScales::Init() {

  static ClassDocumentation<MatchboxTriVecScales> documentation
    ("MatchboxTriVecScales implements scale choices related to transverse momenta and transverse energies for"
     "events with up to three vector bosons in the final state plus additional jets, where the vector bosons"
     "are associated to lepton pairs of distinct families.");

  static Reference<MatchboxTriVecScales,JetFinder> interfaceJetFinder
    ("JetFinder",
     "A reference to the jet finder.",
     &MatchboxTriVecScales::theJetFinder, false, false, true, false, false);

  static Switch<MatchboxTriVecScales,unsigned int> interfaceTriVecScaleChoice
    ("TriVecScaleChoice",
     "The scale choice to use.",
     &MatchboxTriVecScales::theTriVecScaleChoice, 1, false, false);
  static SwitchOption interfaceTriVecScaleChoice1
    (interfaceTriVecScaleChoice,
     "HtPrimeScale",
     "Sum of the transverse energies of the lepton pairs and the transverse momenta of the jets.",
     1);
  static SwitchOption interfaceTriVecScaleChoice2
    (interfaceTriVecScaleChoice,
     "HtPrimeModScale",
     "Sum of the transverse energies of the lepton pairs and the transverse momenta of the jets." 
     "Each transverse jet momentum is thereby suppressed by an exponential of the rapidity difference to the average rapidity of the two hardest jets."
     "This scale choice is defined only for two or more jets in the event.",
     2);
  static SwitchOption interfaceTriVecScaleChoice3
    (interfaceTriVecScaleChoice,
     "EtScale",
     "Sum of the transverse energies of the lepton pairs.",
     3);

}

