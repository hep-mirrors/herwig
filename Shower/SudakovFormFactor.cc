// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "Pythia7/Repository/UseRandom.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/CurrentGenerator.h"

using namespace Herwig;


SudakovFormFactor::~SudakovFormFactor() {}


void SudakovFormFactor::setupLookupTables() {}

void SudakovFormFactor::
get_qz (bool znorm, double p, double R, Energy q0, Energy qmax, Energy &q, double &z) {

  double z0 = .5; 
  Energy qmin; 
  qmin = pow( (pow(q0, 1.+p) - R*pow(100.*GeV, 1.+p))/(1.-R), 1./(1.+p) ); 
  q = pow( pow(qmin, 1.+p) + UseRandom::rnd()
	   *(pow(qmax, 1.+p) - pow(qmin, 1.+p)) , 1./(1.+p) ); 
  z0 = q0/q;

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::get_qz: ==> start extreme <==" << endl;
  }

  if ( q < q0 || z0 >= 0.5) { 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  no branching! " << endl 
			      << "  (qmin < q0 < qmax, q, z0) = ("
			      << qmin << " < " << q0 << " < " << qmax << ", " << q << ", " << z0
			      << ")" << endl;
    }
    q = 0; 
    z = 0;     

  } else {
      if (znorm) {
	// like 1/(1-z)
	z = 1.- (1.-z0)*pow( z0/(1.-z0), UseRandom::rnd() ); 
      } else {
	// like z^2+(1-z)^2 with z0 < z < 1
	do { 
	  z = pow( UseRandom::rnd(), 1./3. );  
	  if ( UseRandom::rndbool() ) z = 1.-z; 
	} while (z < z0 || z > 1. ); 
      }

      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	// generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
	CurrentGenerator::log() << "  branching: (z > z0=q0/q) =  (" 
				<< z << " > " << z0 << ")" 
				<< endl 
				<< "  (qmin < q0 < qmax, q) = ("
				<< qmin << " < " << q0 << " < " << qmax << ", " << q
				<< ")" << endl;
      }
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::get_qz: ==> end extreme <==" << endl;
  }

  return; 
} 
