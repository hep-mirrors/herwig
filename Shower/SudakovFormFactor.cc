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

// inline it!
double SudakovFormFactor::guessz (double z0, double z1) {
  tSplitFun1to2Ptr sF = dynamic_ptr_cast<tSplitFun1to2Ptr>(_splitFun);
  return sF->invIntegOverIntegratedFun( sF->integOverIntegratedFun(z0) + 
     UseRandom::rnd()*( sF->integOverIntegratedFun(z1) - 
			sF->integOverIntegratedFun(z0) ) );
}


Energy2 SudakovFormFactor::guesst (Energy2 t0, Energy2 t1) {
  tSplitFun1to2Ptr sF = dynamic_ptr_cast<tSplitFun1to2Ptr>(_splitFun);
  double z0, z1; 
  z0 = sqrt(t0/t1); 
  z1 = 1.-sqrt(t0/t1); 
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
  Energy2 told, mc2, tmax, t; 
  bool veto = true; 
  tmax = sqr(root_tmax); 
  t = sqr(root_t); 

  //  mc2 = 1.0*GeV2;
  mc2 = sqr(_minScale)/4.;
  t = tmax; 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "SudakovFormFactor::gettz(): extreme ____________________________________________" << endl
			    << "  called with q = " << sqrt(tmax)/GeV 
			    << " and mc = " << sqrt(mc2)/GeV << " (GeV)" << endl; 
  }

  if ( tmax <= 4.*mc2 ) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) 
      CurrentGenerator::log() << "  | tmax < 4mc2! return with (sqrt(t)/GeV, z) = (" 	 
			      << root_t/GeV << ", "
			      << z << ")" << endl; 
    root_t = -1;
    return; 
  }

  // the veto algorithm
  do {
    
    // remind the old value
    told = t; 

    // the larger PS-boundary in z
    z0 = sqrt(mc2/told); 
    z1 = 1.-sqrt(mc2/told); 
    z = guessz(z0, z1); 
    t = guesst(mc2, told); 

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
    z0 = sqrt(mc2/t); 
    z1 = 1.-sqrt(mc2/t); 
    
    // *** ACHTUNG *** collect some sort of statistics of the
    // likelihood of single vetoes in order to check the most likely
    // 1st and only if this fails the 2nd etc 
    veto = false; 
   
    // still inside PS?
    if ((z < z0 || z > z1) && t > 4.*mc2) { 
      veto = true; 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto: not in PS: (z0, z, z1) = ("
				<< z0 << ", "
				<< z << ", "
				<< z1 << ")" << endl;  
      }
      //      psv++; 
      //     cout << "PS (" << vc+1 << ", " << psv << ")" << endl;

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
    if (UseRandom::rnd() > ratio) { 
      veto = true; 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto on P(z)/g(z)" << endl;
      }
      //      hitv++;
      //      cout << "HI (" << vc+1 << ", " << hitv << ")" << endl;
    }

    // alpha_s valid? 
    if ( UseRandom::rnd() > _alpha->value(z*(1.-z)*t)/
	 _alpha->overestimateValue() ) {
      veto = true; 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  X veto on as(q2)/as" << endl;
      }
      // asv++; 
    }

    // alpha_s valid? 
//     if (z*(1.-z)*t < Q2MIN) {
//       veto = true; 
//     };

    // is t valid at all? 
    if (t < 4.*mc2) {
      veto = false; 
      t = -1; 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
	CurrentGenerator::log() << "  | return with no branching, t < 4mc2" << endl;
      }
    }
	
//    if (veto) {
//      vc++; 
//       cout << "(veto, ps, tint, hit, as) = ("
// 	   << vc << ", " 
// 	   << psv << ", " 
// 	   << tintv << ", " 
// 	   << hitv << ", " 
// 	   << asv << ") " 
// 	   << endl; 
//	};
    
//     cout << "veto = " << veto << " (z0, z, z1) = (" 
// 	 << z0 << ", "
// 	 << z << ", "
// 	 << z1 << ")" << endl; 


  } while (veto); 

//   cerr << "---> accepted (t, z) = (" 	 
//        << t << ", "
//        << z << ")" << endl; 
  if (t > 0) root_t = sqrt(t);
  else root_t = -1; ;

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "  return with (sqrt(t)/GeV, z) = (" 	 
			    << root_t/GeV << ", "
			    << z << ")" << endl; 
  }

  return;
}
