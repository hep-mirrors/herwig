// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "Pythia7/Repository/UseRandom.h"

using namespace Herwig;


SudakovFormFactor::~SudakovFormFactor() {}


void SudakovFormFactor::setupLookupTables() {}

void SudakovFormFactor::
get_qz (bool znorm, double p, double R, Energy q0, Energy qmax, Energy &q, double &z) {

  // toy model: returns q according to powerlike (q^p) distribution
  // with cutoff q0, qmin < q0 < qmax.  qmin is chosen such that the
  // probability for a first branching is 1-R.  z is chosen from
  // 1/(1-z) with z0 = m/q and z0 < z < 1-z0; or (if znorm==false)
  // flatter, as z^2+(1-z)^2 with z0 < z < 1.

  double z0 = .5; 
  Energy qmin; 
  qmin = pow( (pow(q0, 1.+p) - R*pow(100.*GeV, 1.+p))/(1.-R), 1./(1.+p) ); 
  q = pow( pow(qmin, 1.+p) + UseRandom::rnd()
	   *(pow(qmax, 1.+p) - pow(qmin, 1.+p)) , 1./(1.+p) ); 

  if (q < q0 ) { 
    q = 0; 
    z = 0; 
  } else {
      z0 = q0/q;
      if (znorm) {
	// like 1/(1-z)
	z = 1. - z0*pow( (1.-z0)/z0, UseRandom::rnd()); 
      } else {
	// like z^2+(1-z)^2 with z0 < z < 1
	do { 
	  z = pow( UseRandom::rnd(), 1./(1.+p) );  
	  if ( UseRandom::rndbool() ) z = 1.-z; 
	} while (z < z0); 
      }
  }
  return; 
} 
