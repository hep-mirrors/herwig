// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FS_QtildaShowerKinematics1to2 class.
//

#include "FS_QtildaShowerKinematics1to2.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/Repository/EventGenerator.h"
#include "Pythia7/Repository/CurrentGenerator.h"

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
//      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 			     << " ===> START DEBUGGING <=== "
// 			     << endl; 
//    }

  // check phase space, ie whether we have a negative radical or not. 
  if ( sqr( dz*dqtilde ) - dm2 < 0 ) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren(): "
			    << "Warning! Cannot reconstruct p_perp, z > sqrt(m2/qtilda2)!" 
			    << endl; 
  }

  // dump all input 
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
			    << "==> start extreme <== " << endl
			    << "  called with (qtilde, z, phi) = (" 
			    << dqtilde << ", " << dz << ", " << dphi << ") " << endl
			    << "  input (alpha, px, py) = (" 
			    << parentSudAlpha << ", " 
			    << parentSudPx << ", " << parentSudPy << ")" << endl; 
  }

  // determine alphas of children according to interpretation of z
  sudAlphaVect.clear(); 
  sudPxVect.clear(); 
  sudPyVect.clear(); 

  sudAlphaVect.push_back( dz*parentSudAlpha ); 
  sudAlphaVect.push_back( (1.-dz)*parentSudAlpha ); 

//   // determine transverse momenta of children
  Energy pPerp; 
  if ( sqr( dz*dqtilde ) - dm2 > 0 ) {
    pPerp = (1.-dz)*sqrt( sqr( dz*dqtilde ) - dm2);
  } else {
    // *** ACHTUNG! *** enforce 'in-phase-space' for the moment..."
    pPerp = 0;
    // *** ACHTUNG! *** replace with proper debug stuff... after
    // proper interfacing.
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren()"
			    << "Warning! Can't get p_perp, " << endl 
			    << " z = "
			    << dz << " < sqrt(m2/qtilda2) = sqrt(" 
			    << dm2 << "/" << sqr(dqtilde) << ") = " 
			    << sqrt(dm2)/dqtilde << "!" << endl;
  }
  sudPxVect.push_back( pPerp*cos(dphi) + dz*parentSudPx ); 
  sudPyVect.push_back( pPerp*sin(dphi) + dz*parentSudPy ); 
  sudPxVect.push_back( - pPerp*cos(dphi) + (1.-dz)*parentSudPx ); 
  sudPyVect.push_back( - pPerp*sin(dphi) + (1.-dz)*parentSudPy ); 

  // dump all output
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
    CurrentGenerator::log()  << "  result (alpha1, alpha2), pPerp = (" 
			     << sudAlphaVect[0] << ", " << sudAlphaVect[1] << "), " 
			     << pPerp << endl
			     << "  result (px1, py1), (px2, py2) = (" 
			     << sudPxVect[0] << ", " << sudPyVect[0] << "), (" 
			     << sudPxVect[1] << ", " << sudPyVect[1] << ")" << endl
			     << "FS_QtildaShowerKinematics1to2::updateChildren() "
			     << "==> end extreme <== " << endl;
  }

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 			    << " ===> END DEBUGGING <=== " << endl; 
//   }

}


void FS_QtildaShowerKinematics1to2::updateChildren( const tShoParPtr theParent, 
						    const CollecShoParPtr theChildren ) {

  if (theChildren.size() != 2) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() " 
			    << "Warning! too many children!" << endl;      
  } else {
    
    double dqtilde = qtilde();
    double dz = z(); 
    double dphi = phi();
    double dm2 = sqr( theParent->data().mass() ); 
    
    // check phase space, ie whether we have a negative radical or not. 
    if ( sqr( dz*dqtilde ) - dm2 < 0 ) {
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren(): "
			      << "Warning! Cannot reconstruct p_perp, z > sqrt(m2/qtilda2)!" 
			      << endl; 
    }
    
    // dump relevant input data
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
			      << "==> start extreme <== " << endl
			      << "  called with (qtilde, z, phi) = (" 
			      << dqtilde << ", " << dz << ", " << dphi << ") " << endl
			      << "  input (alpha, px, py) = (" 
			      << theParent->sudAlpha() << ", " 
			      << theParent->sudPx() << ", " 
			      << theParent->sudPy() << ")" << endl; 
    }

    // determine alphas of children according to interpretation of z
    theChildren[0]->sudAlpha( dz*theParent->sudAlpha() ); 
    theChildren[1]->sudAlpha( (1.-dz)*theParent->sudAlpha() ); 
    
    // determine transverse momenta of children
    Energy pPerp; 
    if ( sqr( dz*dqtilde ) - dm2 > 0 ) {
      pPerp = (1.-dz)*sqrt( sqr( dz*dqtilde ) - dm2);
    } else {
      // *** ACHTUNG! *** enforce 'in-phase-space' for the moment..."
      pPerp = 0;
      // *** ACHTUNG! *** replace with proper debug stuff... after
      // proper interfacing.
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren()"
			      << "Warning! Can't get p_perp, " << endl 
			      << " z = "
			      << dz << " < sqrt(m2/qtilda2) = sqrt(" 
			      << dm2 << "/" << sqr(dqtilde) << ") = " 
			      << sqrt(dm2)/dqtilde << "!" << endl;
    }

    theChildren[0]->sudPx( pPerp*cos(dphi) + dz*theParent->sudPx() );
    theChildren[0]->sudPy( pPerp*sin(dphi) + dz*theParent->sudPy() );
    theChildren[1]->sudPx( - pPerp*cos(dphi) + (1.-dz)*theParent->sudPx() );
    theChildren[1]->sudPy( - pPerp*sin(dphi) + (1.-dz)*theParent->sudPy() );

    // dump all output
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log()  << "  result (alpha1, alpha2), pPerp = (" 
			       << theChildren[0]->sudAlpha() << ", " 
			       << theChildren[1]->sudAlpha() << ", " 
			       << pPerp << endl
			       << "  result (px1, py1), (px2, py2) = (" 
			       << theChildren[0]->sudPx() << ", " 
			       << theChildren[0]->sudPy() << "), (" 
			       << theChildren[1]->sudPx() << ", " 
			       << theChildren[1]->sudPy() << ")" 
			       << endl 
			       << "FS_QtildaShowerKinematics1to2::updateChildren() "
			       << "==> end extreme <== " << endl;
    }
  }
}


void FS_QtildaShowerKinematics1to2::
updateParent( const tShoParPtr theParent, 
	      const CollecShoParPtr theChildren ) {

  if (theChildren.size() != 2) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateParent() " 
			    << "Warning! too many children!" << endl;      
  } else {
    theParent->sudBeta( theChildren[0]->sudBeta() + theChildren[1]->sudBeta() ); 
    theParent->momentum( theChildren[0]->momentum() + theChildren[1]->momentum() ); 

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateParent(): "	 
			      << endl 
			      << "  set beta = " << theParent->sudBeta() 
			      << " = " 
			      << theChildren[0]->sudBeta() 
			      << " + " 
			      << theChildren[1]->sudBeta() << "." 
			      << endl
			      << "  set q = " << theParent->momentum()
			      << endl
			      << "        = " << theChildren[0]->momentum()
			      << endl 
			      << "        + " << theChildren[1]->momentum()
			      << endl
			      << "  q_sud = " 
			      << sudakov2Momentum( theParent->sudAlpha(), 
						   theParent->sudBeta(), 
						   theParent->sudPx(), 
						   theParent->sudPy() )
			      << " from Sudakov variables."
			      << endl; 
    }

  }
  
}


void FS_QtildaShowerKinematics1to2::updateLast( const tShoParPtr theLast ) {

  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  theLast->sudBeta( (sqr(theLast->data().mass()) + theLast->sudPperp2() 
  		     - sqr( theLast->sudAlpha() )*_pVector.m2())
  		    / ( 2.*theLast->sudAlpha()*_pVector*_nVector ) ); 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateLast(): "	 
			    << endl 
			    << "  set beta = " << theLast->sudBeta() 
			    << endl
			    << "  used (q2, alpha, perp2, p.n) = (" 
			    << sqr(theLast->data().mass()) << ", "
			    << theLast->sudAlpha() << ", "
			    << theLast->sudPperp2() << ", "
			    << _pVector*_nVector << ")"
			    << endl; 
  }
      
  // set that new momentum    
  theLast->momentum(  sudakov2Momentum( theLast->sudAlpha(), theLast->sudBeta(), 
					theLast->sudPx(), theLast->sudPy() )  );
}


Energy FS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}



