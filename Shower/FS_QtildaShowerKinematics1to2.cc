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

  // this version is not called! only the next method! 

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
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() extreme ________________________" 
			    << endl
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
			     << sudPxVect[1] << ", " << sudPyVect[1] << ")" << endl;
// 			     << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 			     << "==> end extreme <== " << endl;
  }

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
//     CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() "
// 			    << " ===> END DEBUGGING <=== " << endl; 
//   }

}


void FS_QtildaShowerKinematics1to2::updateChildren( const tShowerParticlePtr theParent, 
						    const ParticleVector theChildren ) {

  if (theChildren.size() != 2) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() " 
			    << "Warning! too many children!" << endl;      
  } else {
    ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
    ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
    double dqtilde = qtilde();
    double dz = z(); 
    double dphi = phi();
    bool glueEmits = false; 
    bool gqq = false;
    if ( theParent->data().id() == 21 ) glueEmits = true; 
    if ( theParent->data().id() == 21 && c1->data().id() == -c2->data().id() ) 
      gqq = true; 

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren() extreme ________________________" << endl
			      << "  found " << theParent->data().PDGName() 
			      << "->" << c1->data().PDGName() 
			      << "+" << c2->data().PDGName()
			      << ", i.e. gqq = " << (gqq ? "true" : "false")
			      << endl
			      << "  called with (qtilde, z, phi) = (" 
			      << dqtilde << ", " << dz << ", " << dphi << ") " << endl
			      << "  input (alpha, px, py) = (" 
			      << theParent->sudAlpha() << ", " 
			      << theParent->sudPx() << ", " 
			      << theParent->sudPy() << ")" << endl; 
    }

    // determine alphas of children according to interpretation of z
    c1->sudAlpha( dz*theParent->sudAlpha() ); 
    c2->sudAlpha( (1.-dz)*theParent->sudAlpha() ); 

    // determine transverse momenta of children
    Energy2 dm2 = Energy2(); 
    Energy pPerp, kinCutoff; 
    // kinCutoff = .15*GeV;
    kinCutoff = kinScale();
    bool inps = true; 
    if (glueEmits) {
      // dm2 = c1->data().mass(); 
      //      dm2 = sqr(max(resScale(), 4.*c1->data().mass()));    
      //      dm2 = sqr(max(kinCutoff, 4.*c1->data().mass()));    
      dm2 = sqr(max(kinCutoff, c1->data().mass()));    
      if ( sqr(dz*(1.-dz)*dqtilde) - dm2 > 0 ) 
	//	pPerp = sqrt(sqr(dz*(1.-dz)*dqtilde) - dm2);
	pPerp = sqrt(sqr(dz*(1.-dz)*dqtilde) - dm2);
      else inps = false;
    } else {
      // dm2 = sqr( theParent->data().mass() ); 
      //      dm2 = sqr(max(resScale(), theParent->data().mass()));
      dm2 = sqr(max(kinCutoff, theParent->data().mass()));
      if ( sqr(1.-dz)*(sqr( dz*dqtilde ) - dm2) > dz*sqr(kinCutoff) ) 
	pPerp = sqrt( sqr(1.-dz)*(sqr( dz*dqtilde ) - dm2) - dz*sqr(kinCutoff) );
      else inps = false; 
    }

    if ( !inps ) {
      // *** ACHTUNG! *** enforce 'in-phase-space' for the moment..."
      pPerp = 0;
      // *** ACHTUNG! *** replace with proper debug stuff... after
      // proper interfacing.
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateChildren()"
			      << " Warning! Can't get p_perp, " << endl 
			      << "  z = "
			      << dz;
      if (glueEmits) {
	CurrentGenerator::log() << ", (m2, qt2) = ("
				<< dm2 << ", " << sqr(dqtilde) 
				<< "), (z(1-z) qt)^2 = " 
				<< sqr(dz*(1.-dz)*dqtilde) << "!" << endl;
	
      } else {
	CurrentGenerator::log() << " < sqrt(m2/qtilda2) = sqrt(" 
				<< dm2 << "/" << sqr(dqtilde) << ") = " 
				<< sqrt(dm2)/dqtilde << "!" << endl;
      }
    }

    c1->sudPx( pPerp*cos(dphi) + dz*theParent->sudPx() );
    c1->sudPy( pPerp*sin(dphi) + dz*theParent->sudPy() );
    c2->sudPx( - pPerp*cos(dphi) + (1.-dz)*theParent->sudPx() );
    c2->sudPy( - pPerp*sin(dphi) + (1.-dz)*theParent->sudPy() );

    // dump all output
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {
      CurrentGenerator::log()  << "  result (alpha1, alpha2), pPerp = (" 
			       << c1->sudAlpha() << ", " 
			       << c2->sudAlpha() << ", " 
			       << pPerp << endl
			       << "  result (px1, py1), (px2, py2) = (" 
			       << c1->sudPx() << ", " 
			       << c1->sudPy() << "), (" 
			       << c2->sudPx() << ", " 
			       << c2->sudPy() << ")" 
			       << endl 
			       << "FS_QtildaShowerKinematics1to2::updateChildren() "
			       << "==> end extreme <== " << endl;
    }
  }
}


void FS_QtildaShowerKinematics1to2::
updateParent( const tShowerParticlePtr theParent, 
	      const ParticleVector theChildren ) {

  if (theChildren.size() != 2) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateParent() " 
			    << "Warning! too many children!" << endl;      
  } else {
    ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
    ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
    theParent->sudBeta( c1->sudBeta() + c2->sudBeta() ); 
    theParent->set5Momentum( c1->momentum() + c2->momentum() ); 

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateParent(): full ____________________________"	 
			      << endl 
			      << "  set beta = " << theParent->sudBeta() 
			      << " = " 
			      << c1->sudBeta() 
			      << " + " 
			      << c2->sudBeta() << "." 
			      << endl
			      << "  set q = " << theParent->momentum()
			      << endl
			      << "        = " << c1->momentum()
			      << endl 
			      << "        + " << c2->momentum()
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


void FS_QtildaShowerKinematics1to2::updateLast( const tShowerParticlePtr theLast ) {

  // set beta component and consequently all missing data from that,
  // using the nominal (i.e. PDT) mass.
  
  Energy theMass; 
  if ( theLast->data().id() == 21 ) theMass = theLast->momentum().mass(); 
  else theMass = theLast->data().mass(); 
  theLast->sudBeta( (sqr(theMass) + theLast->sudPperp2() 
  		     - sqr( theLast->sudAlpha() )*_pVector.m2())
  		    / ( 2.*theLast->sudAlpha()*_pVector*_nVector ) ); 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    CurrentGenerator::log() << "FS_QtildaShowerKinematics1to2::updateLast(): full ______________________________"	 
			    << endl 
			    << "  set beta = " << theLast->sudBeta() 
			    << endl
			    << "  used (sqrt(q2)=m, alpha, perp2, p.n) = (" 
      //			    << sqr(theLast->data().mass()) << ", "
			    << theMass << ", "
			    << theLast->sudAlpha() << ", "
			    << theLast->sudPperp2() << ", "
			    << _pVector*_nVector << ")"
			    << endl; 
  }
      
  // set that new momentum    
  theLast->set5Momentum(  sudakov2Momentum( theLast->sudAlpha(), theLast->sudBeta(), 
					theLast->sudPx(), theLast->sudPy() )  );
}


Energy FS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}



