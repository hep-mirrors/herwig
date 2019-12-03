// -*- C++ -*-
//
// ShowerVariation.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerVariation_H
#define HERWIG_ShowerVariation_H

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {

using namespace ThePEG;

/**
 * A struct identifying a shower variation
 */
struct ShowerVariation {
  
  /**
   * Default constructor
   */
  ShowerVariation() : renormalizationScaleFactor(1.0),
		      factorizationScaleFactor(1.0),
		      firstInteraction(true),
		      secondaryInteractions(false) {}
  
  /**
   * Vary the renormalization scale by the given factor.
   */
  double renormalizationScaleFactor;
  
  /**
   * Vary the factorization scale by the given factor.
   */
  double factorizationScaleFactor;
  
  /**
   * Apply the variation to the first interaction
   */
  bool firstInteraction;
  
  /**
   * Apply the variation to the secondary interactions
   */
  bool secondaryInteractions;
  
  /**
   * Parse from in file command
   */
  string fromInFile(const string&);
  
  /**
   * Put to persistent stream
   */
  void put(PersistentOStream& os) const;
  
  /**
   * Get from persistent stream
   */
  void get(PersistentIStream& is);
  
};

inline PersistentOStream& operator<<(PersistentOStream& os, const ShowerVariation& var) {
  var.put(os); return os;
} 

inline PersistentIStream& operator>>(PersistentIStream& is, ShowerVariation& var) {
  var.get(is); return is;
} 

}
#endif /* HERWIG_ShowerVariation_H */
