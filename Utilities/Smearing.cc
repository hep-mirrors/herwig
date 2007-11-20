// -*- C++ -*-
//
// Smearing.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Smearing class.
//

#include "Smearing.h"
#include <ThePEG/Repository/UseRandom.h>

using namespace Herwig;
// using namespace ThePEG;


bool Smearing::gaussianSmearing(const double mean, 
				const double sigma, 
				double & x ) {
  double xN01, trash; 
  if ( ! azimuthalSmearing( -2.0*log( UseRandom::rnd() ), xN01, trash ) ) 
    return false;
  x = mean + sigma*xN01;
  return true;
}


bool Smearing::azimuthalSmearing(const double rho, 
				 double & vx, 
				 double & vy ) {
  double cosine = 2.0 * UseRandom::rnd() - 1.0;
  double sine   = 2.0 * UseRandom::rnd() - 1.0;
  double cs = sqr(cosine) + sqr(sine);
  if ( cs > 1.0  || cs <= 0.0 ) return false;
  vx = ( sqr(cosine) - sqr(sine) ) * rho/cs;
  vy = 2.0 * cosine * sine * rho/cs;  
  return true;
}

