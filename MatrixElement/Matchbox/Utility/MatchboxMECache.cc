// -*- C++ -*-
//
// MatchboxMECache.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxMECache class.
//

#include "MatchboxMECache.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <numeric>
using std::accumulate;

using namespace Herwig;

MatchboxMECache::MatchboxMECache() {}

MatchboxMECache::~MatchboxMECache() {}

IBPtr MatchboxMECache::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxMECache::fullclone() const {
  return new_ptr(*this);
}

// binary function accumulating hashes
// from Lorentz5Momentum
struct accMomentumHash {

  size_t operator() (size_t before, const Lorentz5Momentum& mom) const {
    // indeed should mask some non-significant digits
    // here as well, but stay with this for the moment
    // as we're having several evaluations with
    // exactly the _same_ ps point 
    size_t res = before;
#ifdef ThePEG_HAS_UNITS_CHECKING
    boost::hash_combine(res,mom.x().rawValue());
    boost::hash_combine(res,mom.y().rawValue());
    boost::hash_combine(res,mom.z().rawValue());
#else
    boost::hash_combine(res,mom.x());
    boost::hash_combine(res,mom.y());
    boost::hash_combine(res,mom.z());
#endif
    return res;
  }

};

size_t MatchboxMECache::hashPhaseSpace() const {
  return accumulate(meMomenta().begin() + 2, meMomenta().end(),
		    0, accMomentumHash());
}

bool MatchboxMECache::calculateME2(double& xme2,
				   const pair<int,int>& corr) {
  map<MECacheKey,double>::const_iterator cached
    = theME2Cache.find(MECacheKey(hashPhaseSpace(),mePartonData(),corr));
  if ( cached == theME2Cache.end() ) {
    xme2 = 0.0;
    return true;
  }
  xme2 = cached->second;
  return false;
}

void MatchboxMECache::cacheME2(double xme2,
			       const pair<int,int>& corr) {
  theME2Cache[MECacheKey(hashPhaseSpace(),mePartonData(),corr)] = xme2;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxMECache::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb;
}

void MatchboxMECache::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxMECache,HandlerBase>
  describeHerwigMatchboxMECache("Herwig::MatchboxMECache", "HwMatchbox.so");

void MatchboxMECache::Init() {

  static ClassDocumentation<MatchboxMECache> documentation
    ("MatchboxMECache provides caching for matrix elements.");

}

