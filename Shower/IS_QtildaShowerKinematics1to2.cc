// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;


IS_QtildaShowerKinematics1to2::~IS_QtildaShowerKinematics1to2() {}


void IS_QtildaShowerKinematics1to2::
updateChildren( const double parentSudAlpha, 
		const Energy parentSudPx, const Energy parentSudPy, 
		vector<double> & sudAlphaVect, 
		vector<Energy> & sudPxVect, vector<Energy> & sudPyVect ) {
  // this is empty...
}


void IS_QtildaShowerKinematics1to2::
updateParent( tCollecShoKinPtr & shoKinChildren ) {
  // this is empty
}

void IS_QtildaShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ParticleVector theChildren ) {
  // this is empty...
}


void IS_QtildaShowerKinematics1to2::
updateParent( const tShowerParticlePtr theParent, 
	      const ParticleVector theChildren ) {

  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
  // get shower variables from 1st child in order to keep notation
  // parent->(c1, c2) clean even though the splitting was initiated
  // from c1.  The name updateParent is still referring to the
  // timelike branching though.

  bool glueEmits = false; 
  if ( theParent->data().id() == 21 ) glueEmits = true; 
  Energy kinCutoff;
  if (glueEmits) {
    kinCutoff = (kinScale() - 0.3*c1->data().mass())/2.3;
  } else {
    kinCutoff = (kinScale() - 0.3*theParent->data().mass())/2.3;
  }
  double mz = 1.-c1->showerKinematics()->z();
  double cphi = cos(c1->showerKinematics()->phi());
  double sphi = sqrt(1.-sqr(cphi));
//   Energy pt = sqrt(sqr((1.-mz)*c1->showerKinematics()->qtilde()) - 
// 		   (1.-mz)*sqr(kinCutoff));
  Energy pt = sqrt(sqr(mz*c1->showerKinematics()->qtilde()) - 
 		   (1.-mz)*sqr(kinCutoff));
  Energy kx = mz*theParent->sudPx() - cphi*pt;
  Energy ky = mz*theParent->sudPy() - sphi*pt; 
  double alpha = theParent->sudAlpha();
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    CurrentGenerator::log()
      << "  updateParent: theParent->sudAlpha() = "
      << theParent->sudAlpha() << endl
      << "  theParent->id() = "
      << theParent->id() 
      << endl
      << "  (z, pt, alpha) = ("
      << 1.-mz << ", "
      << pt << ", " 
      << alpha << ")" << endl
      << "  (cphi, kx, ky) = (" 
      << cphi << ", "
      << kx << ", "
      << ky << ")" << endl;
  }
  // on-shell child
  c2->sudAlpha(mz*alpha);
  c2->sudPx(kx);
  c2->sudPy(ky);
  
  
  // hack to debug IS shower!  Make massive gluons here dependent on
  // OnOffFSShower switch...
  double beta = 0.0;
  if (c2->data().id() == 21) {
    beta = (sqr(750*MeV)+c2->sudPperp2())
      /c2->sudAlpha()/p_dot_n()/2.;
  } else {
    beta = (sqr(c2->data().constituentMass())+c2->sudPperp2())
      /c2->sudAlpha()/p_dot_n()/2.;
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    CurrentGenerator::log() 
      << "  (m, beta, sudPerp2) = (" 
      << c2->data().mass() << ", " 
      << beta << ", " << c2->sudPperp2() << ")" << endl
      << "  (alpha, p.n) = (" << c2->sudAlpha() << ", " 
      << p_dot_n() << ")" << endl
      << "  (m2+pt2, 2 alpha p.n) = (" 
      << sqr(c2->data().mass())+c2->sudPperp2() 
      << ", " << c2->sudAlpha()*p_dot_n()*2. << ")" << endl;
  }
  
  c2->sudBeta(beta);
  c2->set5Momentum(sudakov2Momentum(c2->sudAlpha(), c2->sudBeta(), 
				    c2->sudPx(), c2->sudPy()));
  // spacelike child:
  c1->sudAlpha(theParent->sudAlpha() - c2->sudAlpha());
  c1->sudBeta(theParent->sudBeta() - c2->sudBeta());
  c1->sudPx(theParent->sudPx() - c2->sudPx());
  c1->sudPy(theParent->sudPy() - c2->sudPy());
  c1->set5Momentum(theParent->momentum() - c2->momentum()); 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    CurrentGenerator::log() 
      << "  After 'reconstruction': "
      << "  p->momentum() = " << theParent->momentum() << endl
      << "  c1->momentum() = " << c1->momentum() << endl
      << "  c2->momentum() = " << c2->momentum() 
      << ", m2 = " << c2->momentum().m2() 
      << ", mass = " << c2->momentum().mass() << endl;
  }
}


void IS_QtildaShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast ) {
  //  theLast->sudAlpha(1.0);
  //  cout << "theLast->x() = " << theLast->x() << endl;
  theLast->sudAlpha(theLast->x());
  //  theLast->sudBeta(0.0);
  theLast->sudBeta(sqr(theLast->data().mass())/
		   theLast->sudAlpha()/p_dot_n()/2.);
  theLast->sudPx(0.0*GeV);
  theLast->sudPy(0.0*GeV);
  //  theLast->set5Momentum(_pVector);
  theLast->set5Momentum(
    sudakov2Momentum(theLast->sudAlpha(), theLast->sudBeta(), 
		     theLast->sudPx(), theLast->sudPy()));
}


Energy IS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}
