// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SudakovFormFactor class.
//

#include "SudakovFormFactor.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/Repository/CurrentGenerator.h"

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

// inline it!
double SudakovFormFactor::guessz (double z0, double z1) {
  tSplitFun1to2Ptr sF = dynamic_ptr_cast<tSplitFun1to2Ptr>(_splitFun);
  return sF->invIntegOverIntegratedFun( sF->integOverIntegratedFun(z0) + 
     UseRandom::rnd()*( sF->integOverIntegratedFun(z1) - 
			sF->integOverIntegratedFun(z0) ) );
}


Energy2 SudakovFormFactor::guesst (Energy2 t0, Energy2 t1, double z0, double z1) {
  tSplitFun1to2Ptr sF = dynamic_ptr_cast<tSplitFun1to2Ptr>(_splitFun);
//   double z0, z1; 
//   z0 = sqrt(t0/t1); 
//   z1 = 1.-sqrt(t0/t1); 
//   cerr << "(t1, G(z0), G(z1), as(Q02), exponent) = (" 
//        << t1 << ", " 
//        << sF->integOverIntegratedFun(z1) << ", "
//        << sF->integOverIntegratedFun(z0) << ", "
//        << _alpha->overestimateValue()/(2.*pi) << ", "
//        << 1./( (sF->integOverIntegratedFun(z1) -
// 		sF->integOverIntegratedFun(z0))* 
// 	       _alpha->overestimateValue()/(2.*pi) ) 
//        << ")" << endl; 
  return t1*pow( UseRandom::rnd(), 
		 1./( (sF->integOverIntegratedFun(z1) -
		       sF->integOverIntegratedFun(z0))* 
		      _alpha->overestimateValue()/(2.*pi) )); 
}


void SudakovFormFactor::gettz (Energy root_tmax, Energy &root_t, double &z) {

  tSplitFun1to2Ptr sF = dynamic_ptr_cast<tSplitFun1to2Ptr>(_splitFun);

  double z0, z1, ratio;
  Energy2 told, t0, tmax, t, tmin; 
  bool veto = true; 
  tmax = sqr(root_tmax); 
  t = sqr(root_t); 
//   // Energy kinCutoff = _minScale;
//   Energy kinCutoff = kinScale();
//   bool glueEmits = false; 
//   if (sF->idEmitter() == 21) glueEmits = true; 
//   if ( glueEmits ) { 
//     // have a cutoff of 
//     //    t0 = 16.*sqr(max(_minScale, sF->massFirstProduct()));    
//     //    t0 = sqr(max(_minScale, 4.*sF->massFirstProduct()));    
//     t0 = sqr(max(kinCutoff, sF->massFirstProduct()));    
//   } else {
//     t0 = sqr(max(kinCutoff, sF->massEmitter()));    
//   }

  Energy kinCutoff;
  bool glueEmits = false; 
  if (sF->idEmitter() == 21) glueEmits = true; 
  if ( glueEmits ) { 
    //    kinCutoff = (kinScale() - 0.15*sF->massFirstProduct())/2.3;
    kinCutoff = (kinScale() - 0.3*sF->massFirstProduct())/2.3;
    t0 = sqr(max(kinCutoff, sF->massFirstProduct()));    
  } else {
    //    kinCutoff = (kinScale() - 0.15*sF->massEmitter())/2.3;
    kinCutoff = (kinScale() - 0.3*sF->massEmitter())/2.3;
    t0 = sqr(max(kinCutoff, sF->massEmitter()));    
  }
  t = tmax; 
  tmin = max(t0, sqr(_minScale));

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::gettz(): extreme ____________________________________________" << endl
			    << "  called with q = " << sqrt(tmax)/GeV 
			    << " and q0 = " << sqrt(t0)/GeV << " (GeV)" << endl; 
  }

  if ( tmax <= t0 ) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) 
      CurrentGenerator::log() << "  | tmax < t0! return with (sqrt(t)/GeV, z) = (" 	 
			      << root_t/GeV << ", "
			      << z << ")" << endl; 
    root_t = -1;
    return; 
  }

  //  static long psest, psveto, pgveto, asveto, calls; 
  //  psest = psveto = pgveto = asveto = 0;
  //  calls++;
  // the veto algorithm loop
  do {
    
    // remind the old value
    told = t; 

    // the larger PS-boundary in z (could be part of ShowerKinematics, really!)
    if (glueEmits) {
      z0 = (1.-sqrt(1.-4.*sqrt(t0/told)))/2.;
      z1 = (1.+sqrt(1.-4.*sqrt(t0/told)))/2.;
    } else {
      z0 = sqrt(t0/told)/2.;
      z1 = 1.-kinCutoff/sqrt(told)/2.; // a little overestimate...
    }

    // quick hack gives no more branching without PS
    // if (z0>z1) z1=z0;
    t = guesst(t0, told, z0, z1); 
    z = guessz(z0, z1); 

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "  old->new | (z0, z, z1) = "
			      << sqrt(told)/GeV << "->"
			      << sqrt(t)/GeV << " | "
			      << z0 << ", "
			      << z << ", "
			      << z1 << ")" 
			      << endl; 
    }

    // actual values for z-limits
    if (glueEmits) {
      z0 = (1.-sqrt(1.-4.*sqrt(t0/t)))/2.;
      z1 = (1.+sqrt(1.-4.*sqrt(t0/t)))/2.;
    } else {
      z0 = sqrt(t0/t)/2.;
      z1 = 1.-kinCutoff/sqrt(t)/2.;
      //z1 = 1.-0.1;
    }
    
    // *** ACHTUNG *** collect some sort of statistics of the
    // likelihood of single vetoes in order to check the most likely
    // 1st and only if this fails the 2nd etc 
    veto = false; 
   
    // still inside PS?
    if ((z < z0 || z > z1) && t > tmin) { 
      veto = true; 
      //  psest++;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto: not in PS: (z0, z, z1) = ("
				<< z0 << ", "
				<< z << ", "
				<< z1 << ")" << endl;  
      }
    } else {
      // still REALLY inside PS? (onverestimated allowed PS)
      if (glueEmits) {
	if (sqr(z*(1.-z))*t < t0) {
	  veto = true; 
	  // psveto++;
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	    CurrentGenerator::log() << "  X veto: really not in PS" << endl; 
	  }
	}
      } else {
	if (sqr(1.-z)*(t*sqr(z) - t0) < z*sqr(kinCutoff)) {
	  veto = true; 
	  //  psveto++;
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	    CurrentGenerator::log() << "  X veto: really not in PS" << endl; 
	  }
	}
      } 
    }

    // overestimate of t-integral due to overestimate 
    // of z-integration limits: veto on this!
//     ratio = (int_g_gg(1.-sqrt(MC2/t)) - int_g_gg(sqrt(MC2/t)))
//       /(int_g_gg(1.-sqrt(MC2/told)) - int_g_gg(sqrt(MC2/told))); 
//     if (drand48() > ratio) {
//       veto = true; 
//       tintv++; 
//       //      cout << "TI (" << vc+1 << ", " << tintv << ")" << endl;
//     }; 

    // hit the density? 
    ratio = sF->integratedFun(z, t)/
      sF->overestimateIntegratedFun(z);
    if (!veto && UseRandom::rnd() > ratio) { 
      veto = true; 
      //pgveto++;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto on P(z)/g(z)" << endl;
      }
    }
    // alpha_s valid? 
    Energy2 pt2;
    // simple argument of alpha_s
    pt2 = sqr(z*(1.-z))*t;
    // better argument of alpha_s is the actual pt...
//     if (glueEmits) {
//       pt2 = sqr(z*(1.-z))*t - t0;
//     } else {
//       pt2 = sqr((1.-z))*(sqr(z)*t - t0) - z*sqr(kinCutoff);
//     }

    if ( !veto && UseRandom::rnd() > _alpha->value(pt2)/
	 _alpha->overestimateValue() ) {
      veto = true; 
      //asveto++;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto on as(q2)/as = " 
				<< _alpha->value(pt2)
				<< "/" << _alpha->overestimateValue() 
				<< " = " 
				<< _alpha->value(pt2)/_alpha->overestimateValue() 
				<< endl;
      }
    }
    // is t valid at all? 
    if (t < tmin) {
      //    if (t < t0) {
      veto = false; 
      t = -1; 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  | return with no branching, t < t0" << endl;
      }
    }
  } while (veto); 

//   cerr << "---> accepted (t, z) = (" 	 
//        << t << ", "
//        << z << ")" 
//  cerr << "(calls, ps-est, ps, P/g, as, sum) \\ (per call) = " 
//       << endl << flush; 
//   cerr << "(" << calls << ", "
//        << psest << ", "
//        << psveto << ", "
//        << pgveto << ", "
//        << asveto << ", " 
//        << psest+psveto+pgveto+asveto << ") ("
//        << double(psest)/double(calls) << ", "
//        << double(psveto)/double(calls) << ", "
//        << double(pgveto)/double(calls) << ", "
//        << double(asveto)/double(calls) << ", " 
//        << double(psest+psveto+pgveto+asveto)/double(calls) << ")\r";
// ran this for a good number of events. Variations should only appear with
// crazy parameters...

// (4211, 268, 2271, 1459, 1147, 5145) (0.0636428, 0.539302, 0.346474, 0.272382, 1.2218)


  if (t > 0) root_t = sqrt(t);
  else root_t = -1; ;

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "  return with (sqrt(t)/GeV, z) = (" 	 
			    << root_t/GeV << ", "
			    << z << ")" << endl; 
  }

  return;
}
