// -*- C++ -*-
//
// SU2Helper.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SU2Helper_H
#define HERWIG_SU2Helper_H

#include "ThePEG/PDT/ParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Ken Arnold, Simon Platzer
 *
 * \brief Helpers for book keeping in electroweak processes.
 *
 * @see \ref MatchboxMEBaseInterfaces "The interfaces"
 * defined for MatchboxMEBase.
 */
struct SU2Helper {

  /**
   * Return true, if the left-(right-)handed projection of 
   * this particle (antiparticle) belongs to a weak 
   * SU(2) doublet.
   */
  static bool isInSU2Doublet(tcPDPtr p) {
    return abs(p->id()) < 9 || (abs(p->id()) >= 11 && abs(p->id()) < 19);
  }

  /**
   * Return true, if the left-(right-)handed projection of 
   * this particle (antiparticle) is the up component of a weak 
   * SU(2) doublet.
   */
  static bool isSU2Up(tcPDPtr p) { return abs(p->id())%2==0 && isInSU2Doublet(p); }

  /**
   * Return true, if the left-(right-)handed projection of 
   * this particle (antiparticle) is the down component of a weak 
   * SU(2) doublet.
   */
  static bool isSU2Down(tcPDPtr p) { return abs(p->id())%2!=0 && isInSU2Doublet(p); }

  /**
   * Return the conjugate component in the same weak SU(2) 
   * doublet, or a null pointer if the left-(right-)handed
   * projection of this particle (antiparticle) does not belong 
   * to a weak SU(2) doublet.
   */
  static tcPDPtr SU2CC(tcPDPtr p, int familyShift = 0);

  /**
   * Return the family
   */
  static int family(tcPDPtr p);

};


}

#endif // HERWIG_SU2Helper_H
