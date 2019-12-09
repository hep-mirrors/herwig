// -*- C++ -*-
//
// stat.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Stat_H
#define HERWIG_Stat_H

namespace Herwig{

/**
 * Documentation for the statistics struct. 
 * Copied from ThePEG::StandardEventHandler.
 */
struct Stat {

  /** Standard constructor */
  Stat() : attempted(0), accepted(0), sumw(0.0), maxXSec(CrossSection()),
       totsum(0.0) {}

  /** Constructor with arguments.*/
  Stat(long att, long acc, double w, CrossSection x, double sumw)
    : attempted(att), accepted(acc), sumw(w), maxXSec(x),
       totsum(sumw) {}

  /**
   * Calculation of the cross section.
   */
  inline CrossSection xSec() const {
    return totsum >0.0? maxXSec*sumw/totsum: maxXSec;
  }

  /** Store the number of attempts */
  long attempted;
  /** Store the number of accepted ones */
  long accepted;
  /** Sum of weights */
  double sumw;
  /** Maximal cross section */
  CrossSection maxXSec;
  /** Maximal weights */
  double totsum;

  /** Overloaded += operator */
  const Stat & operator+= (const Stat & s) {
    attempted += s.attempted;
    accepted += s.accepted;
    sumw += s.sumw;
    totsum = max(totsum, s.totsum);
    if ( totsum > 0.0 )
      maxXSec = max(maxXSec, s.maxXSec);
    else
      maxXSec += s.maxXSec;
    return *this;
  }
};
}
#endif
