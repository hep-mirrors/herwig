// -*- C++ -*-
//
// MatchboxTriVecScales.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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

Energy2 MatchboxTriVecScales::HtPrimeScale() const {

  tcPDVector pd (mePartonData().begin() + 2, mePartonData().end());
  vector<LorentzMomentum> p (meMomenta().begin() + 2, meMomenta().end());
  tcPDPtr t1 = mePartonData()[0];
  tcPDPtr t2 = mePartonData()[1];
  tcCutsPtr cuts = lastCutsPtr();

  theJetFinder->cluster(pd, p, cuts, t1, t2);

  Energy sumpartpt = ZERO;
  int foundpart = 0;
  Energy vec1et = ZERO;
  Energy vec2et = ZERO;
  Energy vec3et = ZERO;
  Energy sumvecet = ZERO;
  int foundlept[3] = {0,0,0}; // First entry for e-like, second entry for mu-like, third entry for tau-like
  int ivec1 = 0;
  int ivec2 = 0;
  int ivec3 = 0;
  tcPDVector::const_iterator itpd = pd.begin();
  int ip = 2;
  for (vector<LorentzMomentum>::const_iterator itp = p.begin() ;
       itp != p.end(); ++itp, ++itpd, ++ip ) {
    if ( (**itpd).coloured() ) {
      sumpartpt += (*itp).perp();
      foundpart++;
    }
    if ( !(**itpd).coloured() ) {
      if ( abs((**itpd).id())==ParticleID::nu_e || abs((**itpd).id())==ParticleID::eminus ) {
        if ( foundlept[0]==0 ) {
          foundlept[0] +=1;
          ivec1 = ip;
	}
        if ( foundlept[0]==1 ) {
          vec1et = sqrt( (p[ivec1]+p[ip]).m2() + (p[ivec1]+p[ip]).perp2() );
	}
      }
      if ( abs((**itpd).id())==ParticleID::nu_mu || abs((**itpd).id())==ParticleID::muminus ) {
        if ( foundlept[1]==0 ) {
          foundlept[1] +=1;
          ivec2 = ip;
	}
        if ( foundlept[1]==1 ) {
          vec2et = sqrt( (p[ivec2]+p[ip]).m2() + (p[ivec2]+p[ip]).perp2() );
	}
      }
      if ( abs((**itpd).id())==ParticleID::nu_tau || abs((**itpd).id())==ParticleID::tauminus ) {
        if ( foundlept[2]==0 ) {
          foundlept[2] +=1;
          ivec3 = ip;
	}
        if ( foundlept[2]==1 ) {
          vec3et = sqrt( (p[ivec3]+p[ip]).m2() + (p[ivec3]+p[ip]).perp2() );
	}
      }
    }
  }
  sumvecet = vec1et+vec2et+vec3et;

  // Check for consistency in number of lepton pairs and members therein
  for (int i = 0; i<3; ++i ) {
    if (foundlept[i]>2 || foundlept[i]%2!=0) 
      throw Exception() << "MatchboxTriVecScales::HtPrimeScale(): Inconsistency in number of lepton pairs and members therein!"
			<< Exception::runerror;
  }

  return sqr(sumpartpt+sumvecet);

}

Energy2 MatchboxTriVecScales::HtPrimeModScale() const {
  throw Exception() << "MatchboxTriVecScales::HtPrimeModScale(): Not yet implemented!" << Exception::runerror;
}

Energy2 MatchboxTriVecScales::EtScale() const {
  throw Exception() << "MatchboxTriVecScales::EtScale(): Not yet implemented!" << Exception::runerror;
}

Energy2 MatchboxTriVecScales::renormalizationScale() const {
//   if ( theTriVecScaleChoice==1 ) return 0.5*HtPrimeScale();
//   else if ( theTriVecScaleChoice==2 ) return 0.5*HtPrimeModScale();
//   else if ( theTriVecScaleChoice==3 ) return 0.5*EtScale();
  if ( theTriVecScaleChoice==2 ) return 0.5*HtPrimeModScale();
  if ( theTriVecScaleChoice==3 ) return 0.5*EtScale();
  return 0.5*HtPrimeScale(); // The default is theTriVecScaleChoice==1
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
     "events with up to three vector bosons in the final state, plus additional jets, where the vector bosons"
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
     "Each transverse jet momenta is thereby suppressed by an exponential of the rapidity difference to the average rapidity of the two hardest jets."
     "This scale choice is of benefit only for two or more jets in the event.",
     2);
  static SwitchOption interfaceTriVecScaleChoice3
    (interfaceTriVecScaleChoice,
     "EtScale",
     "Sum of the transverse energies of the lepton pairs and the transverse energies of the jets."
     "This scale choice is of benefit only for two or more jets in the event.",
     3);

}

