// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Smearing class.
//

#include "Smearing.h"
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;
// using namespace Pythia7;


bool Smearing::gaussianSmearing( const double mean, const double sigma, double & x ) {
  double xN01, trash; 
  if ( ! azimuthalSmearing( -2.0*log( UseRandom::rnd() ), xN01, trash ) ) return false;
  x = mean + sigma*xN01;
  return true;
}


bool Smearing::azimuthalSmearing( const double rho, double & vx, double & vy ) {
  double cosine = 2.0 * UseRandom::rnd() - 1.0;
  double sine   = 2.0 * UseRandom::rnd() - 1.0;
  double cs = pow(cosine,2) + pow(sine,2);
  if ( cs > 1.0  || cs <= 0.0 ) return false;
  vx = ( pow(cosine,2) - pow(sine,2) ) * rho/cs;
  vy = 2.0 * cosine * sine * rho/cs;  
  return true;
}

