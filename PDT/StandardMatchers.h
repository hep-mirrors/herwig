// -*- C++ -*-
//
// StandardMatchers.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_StandardMatchers_H
#define Herwig_StandardMatchers_H
// This is the declaration of the AnyMatcher,


#include "ThePEG/PDT/Matcher.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace Herwig {
using namespace ThePEG;

/** \file Herwig++/PDT/StandardMatchers.h
 *
 * This file declare a set of standard matcher classes in addition to those
 * defined in ThePEG. The classes can be used by themselves (with
 * their static functions) or together with the Matcher class to
 * define Interfaced objects of the MatcherBase type to be used in the
 * Repository. Suitable typedefs are declared for the latter.
 *
 * @see Matcher
 * @see MatcherBase
 */

/**
 * A Matcher class which matches photons
 */
struct PhotonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef PhotonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return pd.id()==ParticleID::gamma;
  }
  /** A simplified but unique class name. */
  static string className() { return "Photon"; }
};

/**
 * Gives a MatcherBase class based on PhotonMatcher. 
 */
typedef Matcher<PhotonMatcher> MatchPhoton;

/**
 * A Matcher class which matches top quarks
 */
struct TopMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef TopMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return abs(pd.id())==ParticleID::t;
  }
  /** A simplified but unique class name. */
  static string className() { return "Top"; }
};

/**
 * Gives a MatcherBase class based on TopMatcher. 
 */
typedef Matcher<TopMatcher> MatchTop;

/**
 * A Matcher class which matches any hadron.
 */
struct HadronMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef HadronMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    if(pd.id()!=ParticleID::gamma) return Check(pd.id()); 
    else {
      Ptr<BeamParticleData>::const_pointer beam = 
	dynamic_ptr_cast< Ptr<BeamParticleData>::const_pointer>(&pd);
      if(!beam) return false;
      if(!beam->pdf()) return false;
      return true;
    }
  }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    bool hadron = 
      ((id/10)%10 && (id/100)%10 && (id/1000)%10 == 0) ||
      ((id/10)%10 && (id/100)%10 && (id/1000)%10);
    if(hadron) return true;
    // special for gamma when acting like a hadron
    if(id!=ParticleID::gamma) return false;
    tcPDPtr gamma = CurrentGenerator::current().getParticleData(ParticleID::gamma);
    Ptr<BeamParticleData>::const_pointer beam = 
      dynamic_ptr_cast< Ptr<BeamParticleData>::const_pointer>(gamma);
    if(!beam) return false;
    if(!beam->pdf()) return false;
    return true;
  }
  /** A simplified but unique class name. */
  static string className() { return "Hadron"; }
};
/** Gives a MatcherBase class based on HadronMatcher. */
typedef Matcher<HadronMatcher> MatchHadron;

/**
 * A Matcher class which matches W bosons
 */
struct WBosonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef WBosonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return abs(pd.id())==ParticleID::Wplus;
  }
  /** A simplified but unique class name. */
  static string className() { return "WBoson"; }
};

/**
 * Gives a MatcherBase class based on WBosonMatcher. 
 */
typedef Matcher<WBosonMatcher> MatchWBoson;

/**
 * A Matcher class which matches Z bosons
 */
struct ZBosonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef ZBosonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return abs(pd.id())==ParticleID::Z0;
  }
  /** A simplified but unique class name. */
  static string className() { return "ZBoson"; }
};

/**
 * Gives a MatcherBase class based on ZBosonMatcher. 
 */
typedef Matcher<ZBosonMatcher> MatchZBoson;

/**
 * A Matcher class which matches Higgs bosons
 */
struct HiggsBosonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef HiggsBosonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return abs(pd.id())==ParticleID::h0;
  }
  /** A simplified but unique class name. */
  static string className() { return "HiggsBoson"; }
};

/**
 * Gives a MatcherBase class based on HiggsBosonMatcher. 
 */
typedef Matcher<HiggsBosonMatcher> MatchHiggsBoson;
}

#endif /* Herwig_StandardMatchers_H */
