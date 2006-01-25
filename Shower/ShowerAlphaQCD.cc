// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ShowerAlphaQCD class.
//

#include "ShowerAlphaQCD.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ShowerAlphaQCD.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ShowerAlphaQCD::~ShowerAlphaQCD() {}

void ShowerAlphaQCD::persistentOutput(PersistentOStream & os) const {
  os << _asType << _Qmin;
}

void ShowerAlphaQCD::persistentInput(PersistentIStream & is, int) {
  is >> _asType >> _Qmin;
}

ClassDescription<ShowerAlphaQCD> ShowerAlphaQCD::initShowerAlphaQCD;
// Definition of the static class description member.

void ShowerAlphaQCD::Init() {

  static ClassDocumentation<ShowerAlphaQCD> documentation
    ("This (concrete) class describes the QCD alpha running.");

  static Switch<ShowerAlphaQCD, int> intAsType
    ("NPAlphaS",
     "Behaviour of AlphaS in the NP region",
     &ShowerAlphaQCD::_asType, 1, false, false);
  static SwitchOption intAsTypeZero
    (intAsType, "AsTypeZero","zero below Q_min", 1);
  static SwitchOption intAsTypeConst
    (intAsType, "AsTypeConst","const as(Qmin) below Q_min", 2);
  static SwitchOption intAsTypeLin
    (intAsType, "AsTypeLin ","growing linearly below Q_min", 3);
  static SwitchOption intAsTypeQuad
    (intAsType, "AsTypeQuad","growing quadratically below Q_min", 4);
  static SwitchOption intAsTypeExx1
    (intAsType, "AsTypeExx1 ", "quad from 100 down to as(Q_min)", 5);
  static SwitchOption intAsTypeExx2
    (intAsType, "AsTypeExx2 ", "const = 100 below Q_min", 6);

  // default such that as(Qmin) = 1 in the current parametrization.
  // min = Lambda3
  static Parameter<ShowerAlphaQCD,Energy> intQmin
    ("Qmin", "Q < Qmin is treated with NP parametrization of as (unit [GeV])",
     &ShowerAlphaQCD::_Qmin, GeV, 0.630882*GeV, 0.330445*GeV,
     100.0*GeV,false,false,false);
}

double ShowerAlphaQCD::value(const Energy2 scale) {
  //  OR CALL THE ONE DEFINED IN PYTHIA7 (see class ThePEG::O1AlphaS)
  //  double val = alpha_s(scale, sqr(0.33197*GeV), _asType); gives xe+15...
  // chosen alpha_s(.39GeV) = 176.3 as a maximum value...
  //  double val = alpha_s(scale, sqr(0.630882*GeV), _asType); // gives as = 1
  double val = alpha_s(scale, sqr(_Qmin), _asType); // gives as = 1
  return scaleFactor() * val;
}

double ShowerAlphaQCD::overestimateValue() {
  double val = 0.0; 
  //  if ( _asType < 5 ) val = value(sqr(0.33197*GeV)); 
  //  if ( _asType < 5 ) val = value(sqr(0.5*GeV)); 
  if ( _asType < 5 ) val = value(sqr(_Qmin)); // gives as = 1
  else val = 100.; 
  return scaleFactor() * val; 

}

//////////////////////////////////////////////////////////////////////////
// private stuff:

double ShowerAlphaQCD::alphaTwoLoop(Energy q, Energy lam, short nf) {
  double x, b0, b1, b2; 
  x = sqr(q/lam);
  b0 = 11. - 2./3.*nf;
  b1 = 51. - 19./3.*nf;
  b2 = 2857. - 5033./9.*nf + 325./27.*sqr(nf);
  return( 4.*pi/(b0*log(x))*
	  (1. - 2.*b1/sqr(b0)*log(log(x))/log(x) + 
	   4.*sqr(b1)/(sqr(sqr(b0))*sqr(log(x)))*
	   (sqr(log(log(x)) - 0.5) + b2*b0/(8.*sqr(b1)) - 5./4.)) );  
}

pair<short, Energy> ShowerAlphaQCD::getLamNfTwoLoop(Energy q) {

  // hacked in masses by hand for the moment before proper
  // interfacing...  obtained lambda solutions numerically in
  // Mathematica with my alphas.m using two-loop alphas from PDT2002
  // and as(M_Z=91.187GeV) = 0.118 *** ACHTUNG! *** this HAS to be done
  // automatically acc to the masses and as(M_Z) given by the PDT
  // class (which is supposed to be up-to-date)

  Energy mt, mb, mc;
  mt = 175.0*GeV;
  mb = 4.5*GeV;
  mc = 1.35*GeV;
  Energy lambda3, lambda4, lambda5, lambda6;
  lambda3 = 0.330445*GeV;
  lambda4 = 0.289597*GeV;
  lambda5 = 0.208364*GeV; 
  lambda6 = 0.0880617*GeV;
  short nf = 3;
  Energy lambda = 0.1*GeV;
  
  // get lambda and nf according to the thresholds
  if(q < mc) {
    lambda = lambda3;
    nf = 3;
  } else {
    if(q < mb) {
      lambda = lambda4;
      nf = 4;
    } else {
      if(q < mt) {
	lambda = lambda5;
	nf = 5;
      } else {
	lambda = lambda6;
	nf = 6;
      }
    }
  }  
  return pair<short,Energy>(nf, lambda);
}


double ShowerAlphaQCD::alpha_s(Energy2 q2, Energy2 q2min, int type) {

  pair<short,Energy> nflam;
  Energy q, qmin; 
  q = sqrt(q2); 
  qmin = sqrt(q2min);
  double val = 0.0; 

  nflam = getLamNfTwoLoop(1.*GeV);  // somewhere btw m_s and m_c gives Lambda3
  if ( qmin < nflam.second ) {
    cerr << "alpha_s:  WARNING! q2min chosen smaller than lambda3! return -1" << endl;
    return(-1.); 
  }

  if (q < qmin) {
    nflam = getLamNfTwoLoop(qmin); 
    double val0 = alphaTwoLoop(qmin, nflam.second, nflam.first);
    switch (type) {
    case 1: 
      // flat, zero; the default type with no NP effects.
      val = 0.; 
      break; 
    case 2: 
      // flat, non-zero alpha_s = alpha_s(q2min).
      val = val0;
      break; 
    case 3: 
      // linear in q
      val = val0*q/qmin;
      break; 
    case 4:
      // quadratic in q
      val = val0*q2/q2min;
      break; 
    case 5:
      // quadratic in q, starting off at 100, ending on as(qmin)
      val = (val0 - 100.)*q2/q2min + 100.;
      break; 
    case 6:
      // just big and constant
      val = 100.*val0;
      break; 
    }
  } else {
    // the 'ordinary' case    
    nflam = getLamNfTwoLoop(q); 
    val = alphaTwoLoop(q, nflam.second, nflam.first);
  }; 
  
  return (val); 
}



// ***ACHTUNG*** put somthing more serious here, later
// a parametrization of alpha_s for test purposes only

// double ShowerAlphaQCD::alpha_s(Energy2 q2, Energy2 q2min, int type) {

  // hier mit 4 flavours
  //  Energy lambda = 1.0;
  //  double val = 0.0; 
  //  Energy lambda4 = 0.3;
  // int cflavour;

  // assign different lambda_QCD values for different numbers of 
  // active flavours in a smooth way

  // my old version of this...
//   if(q2 < 1.0*GeV2) {
//     lambda = 0.3280*GeV;
//     cflavour = 27;
//   } else {
//     if(q2 < sqr(5.0*GeV)) {
//       lambda = lambda4;
//       cflavour = 25;
//     } else {
//       if(q2 < sqr(100.*GeV)) {
// 	lambda = 0.2349*GeV;
// 	cflavour = 23;
//       } else {
// 	lambda = 0.1320*GeV;
// 	cflavour = 21;
//       }
//     }
//   };
  
//   // different choices for q2 in the non-perturbative region, q2 < q2min

//   if (q2 < q2min) {
//     if (q2 < 0) {
//       cerr << "alpha_s: negative q2" << endl; 
//       val = -1.;
//     } else { 
//       switch (type) {
//       case 1: 
// 	// flat, zero; the default type with no NP effects.
// 	val = 0.; 
// 	break; 
//       case 2: 
// 	// flat, non-zero alpha_s = alpha_s(q2min).
// 	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)));
// 	break; 
//       case 3: 
// 	// linear 
// 	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)))
// 	  *q2/q2min;
// 	break; 
//       case 4:
// 	// quadratic
// 	val = 12.*M_PI/(((double) cflavour)*log(q2min/sqr(lambda)))
// 	  *sqr(q2/q2min);
// 	break; 
//       }; 
//     }
//   } else {
//     // the common running coupling
//     val = 12.*M_PI/(((double) cflavour)*log(q2/sqr(lambda)));
//   }; 
// return(val); 

// }
