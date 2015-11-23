// -*- C++ -*-
//
// MultiplicityInfo.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respe
#include "Herwig/Utilities/Statistic.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Analysis
 *  Enumeration for species of particle
 */
enum ParticleSpecies {
  lightMeson=0,strangeMeson,lightBaryon,other
};

/** \ingroup Analysis
 *  Struct for the multiplcity data
 */
struct MultiplicityInfo {

  /**
   *  Default constructor
   * @param mult  The observed multiplcity.
   * @param error The error on the observed multiplicity
   * @param type  The type of particle
   */
  MultiplicityInfo(double mult=0.,double error=0.,
		   ParticleSpecies type=other)
    : obsMultiplicity(mult), obsError(error), type(type) {};

  /**
   *  The observed multiplicity
   */
  double obsMultiplicity;

  /**
   *  The error on the observed multiplicity
   */
  double obsError;

  /**
   *  The type of particle
   */
  ParticleSpecies type;

  /**
   *  Simulation statistics for particles of this type
   */
  Statistic count;

  /**
   *  The average number per event
   */
  double simMultiplicity() { return count.mean();}

  /**
   *  The error on the average number per event
   */
  double simError() { return count.mean_stdDev();}

  /**
   *  Is the result more than \f$3\sigma\f$ from the experimental result
   */
  double nSigma() {
    return obsMultiplicity == 0.0 ? 0.0 :
      (simMultiplicity() - obsMultiplicity) 
      / sqrt(sqr(simError()) + sqr(obsError));
  }

  /**
   * Plot standard error in a simple barchart
   */
  string bargraph() {
    if (obsMultiplicity == 0.0) return "     ?     ";
    else if (nSigma() >= 6.0)   return "-----|---->";
    else if (nSigma() >= 5.0)   return "-----|----*";
    else if (nSigma() >= 4.0)   return "-----|---*-";
    else if (nSigma() >= 3.0)   return "-----|--*--";
    else if (nSigma() >= 2.0)   return "-----|-*---";
    else if (nSigma() >= 1.0)   return "-----|*----";
    else if (nSigma() > -1.0)   return "-----*-----";
    else if (nSigma() > -2.0)   return "----*|-----";
    else if (nSigma() > -3.0)   return "---*-|-----";
    else if (nSigma() > -4.0)   return "--*--|-----";
    else if (nSigma() > -5.0)   return "-*---|-----";
    else if (nSigma() > -6.0)   return "*----|-----";
    else                        return "<----|-----";
  }
};
}
