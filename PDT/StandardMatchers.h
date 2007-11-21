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
}
#endif /* Herwig_StandardMatchers_H */
