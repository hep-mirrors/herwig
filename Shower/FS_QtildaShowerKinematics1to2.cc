// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtildaShowerKinematics1to2 class.
//

#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"

using namespace Herwig;


FS_QtildaShowerKinematics1to2::~FS_QtildaShowerKinematics1to2() {}


void FS_QtildaShowerKinematics1to2::
updateChildren( const double parentSudAlpha, const Energy parentSudPx, 
		const Energy parentSudPy, vector<double> & sudAlphaVect,
		vector<Energy> & sudPxVect, vector<Energy> & sudPyVect ) {

  double dqtilde = qtilde();
  double dz = z(); 
  double dphi = phi();
  double dm2 = pVector().mass2(); 

//    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//      generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
//  		       << " ===> START DEBUGGING <=== "
//        //		       << "   EventNumber=" << generator()->currentEventNumber() 
//  		       << endl; 
//      // check phase space, ie whether we have a negative radical or not. 
//      if ( sqr( dz*dqtilde ) - sqr(dm2) < 0 ) {
//        generator()->log() << "Warning! Cannot reconstruct p_perp, z > sqrt(m2/qtilda2)!" << endl; 
//      }
//    }

//   // dump all input 
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
//     generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 		       << " ===> START EXTREME DEBUGGING <=== " << endl
// 		       << "Called with (qtilde, z, phi) = (" 
// 		       << dqtilde << ", " << dz << ", " << dphi << ") " << endl
// 		       << "Input (alpha, px, py) = (" 
// 		       << parentSudAlpha << ", " 
// 		       << parentSudPx << ", " << parentSudPy << ")" << endl; 
//   }

  // determine alphas of children according to interpretation of z
//   sudAlphaVect.clear(); 
//   sudPxVect.clear(); 
//   sudPyVect.clear(); 

//   sudAlphaVect.push_back( dz*parentSudAlpha ); 
//   sudAlphaVect.push_back( (1.-dz)*parentSudAlpha ); 

//   // determine transverse momenta of children
//   Energy pPerp; 
//   if ( sqr( dz*dqtilde ) - sqr(dm2) > 0 ) {
//     pPerp = (1.-dz)*sqrt( sqr( dz*dqtilde ) - sqr(dm2));
//   } else {
//     // *** ACHTUNG! *** enforce 'in-phase-space' for the moment..."
//     pPerp = 0;
//     // *** ACHTUNG! *** replace with proper debug stuff... after
//     // proper interfacing.
//     cerr << "S_QtildaShowerKinematics1to2::updateChildren()"
// 	 << "Warning! Cannot reconstruct p_perp, z > sqrt(m2/qtilda2)!" << endl;
//   }
//   sudPxVect.push_back( pPerp*cos(dphi) + dz*parentSudPx ); 
//   sudPyVect.push_back( pPerp*sin(dphi) + dz*parentSudPy ); 
//   sudPxVect.push_back( - pPerp*cos(dphi) + (1.-dz)*parentSudPx ); 
//   sudPyVect.push_back( - pPerp*sin(dphi) + (1.-dz)*parentSudPy ); 

  // dump all output
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
//     generator()->log() << "result pPerp = " 
// 		       << pPerp << endl
// 		       << "result (alpha1, alpha2) = (" 
// 		       << sudAlphaVect[0] << ", " << sudAlphaVect[1] << ")" << endl
// 		       << "result (px1, py1), (px2, py2) = (" 
// 		       << sudPxVect[0] << ", " << sudPyVect[0] << "), (" 
// 		       << sudPxVect[1] << ", " << sudPyVect[1] << ")" << endl
// 		       << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 		       << " ===> END EXTREME DEBUGGING <=== " << endl;
//   }

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     generator()->log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 		       << " ===> END DEBUGGING <=== " << endl; 
//   }

}


void FS_QtildaShowerKinematics1to2::
updateParent( tCollecShoKinPtr & shoKinChildren ) {

  //***LOOKHERE*** WRITE THE CODE

}


Energy FS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}


