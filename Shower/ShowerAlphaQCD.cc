// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQCD class.
//

#include "ShowerAlphaQCD.h"
#include "Pythia7/Interface/ClassDocumentation.h"

using namespace Herwig;


ShowerAlphaQCD::~ShowerAlphaQCD() {}


ClassDescription<ShowerAlphaQCD> ShowerAlphaQCD::initShowerAlphaQCD;
// Definition of the static class description member.


void ShowerAlphaQCD::Init() {

  static ClassDocumentation<ShowerAlphaQCD> documentation
    ("This (concrete) class describes the QCD alpha running.");

}


double ShowerAlphaQCD::value(const Energy2 scale) {

  //  double val = alpha_s(scale, 1.0*GeV2, 1);
  // q2min could as well be independent of that scale! 
  double val = alpha_s(scale, sqr(_pointerShowerConstrainer->cutoffQScale(ShowerIndex::QCD)), 1);

  // ***ACHTUNG*** just a call to the simpler dummy function used
  // previously in SG's fragmentation programs

  // ***LOOKHERE*** HERE THE FUNCTION alphaQCD( scale ) SHOULD BE DEFINED
  //                OR CALL THE ONE DEFINED IN PYTHIA7
  //                (see class Pythia7::O1AlphaS)

  return scaleFactor() * val;
}


//////////////////////////////////////////////////////////////////////////
// private stuff:


// ***ACHTUNG*** put somthing more serious here, later
// a parametrization of alpha_s for test purposes only

double ShowerAlphaQCD::alpha_s(Energy2 q2, Energy2 q2min, int type) {

  // hier mit 4 flavours
  Energy lambda = 1.0;
  double val = 0.0; 
  Energy lambda4 = 0.3;
  int cflavour;

  // assign different lambda_QCD values for different numbers of 
  // active flavours in a smooth way

  if(q2 < 1.0*GeV2) {
    lambda = 0.3280*GeV;
    cflavour = 27;
  } else {
    if(q2 < sqr(5.0*GeV)) {
      lambda = lambda4;
      cflavour = 25;
    } else {
      if(q2 < sqr(100.*GeV)) {
	lambda = 0.2349*GeV;
	cflavour = 23;
      } else {
	lambda = 0.1320*GeV;
	cflavour = 21;
      }
    }
  };

  // different choices for q2 in the non-perturbative region, q2 < q2min

  if (q2 < q2min) {
    if (q2 < 0) {
      cerr << "alpha_s: negative q2" << endl; 
      val = -1.;
    } else { 
      switch (type) {
      case 1: 
	// flat, zero; the default type with no NP effects.
	val = 0.; 
	break; 
      case 2: 
	// flat, non-zero alpha_s = alpha_s(q2min).
	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)));
	break; 
      case 3: 
	// linear 
	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)))
	  *q2/q2min;
	break; 
      case 4:
	// quadratic
	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)))
	  *sqr(q2/q2min);
	break; 
      }; 
    }
  } else {
    // the common running coupling
    val = 12.*M_PI/(((double) cflavour)*log(q2/sqr(lambda)));
  }; 
  
  return (val); 
}

