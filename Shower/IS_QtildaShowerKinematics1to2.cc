// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"

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
  Energy pt = sqrt(sqr(mz*c1->showerKinematics()->qtilde()) - 
		   (1.-mz)*sqr(kinCutoff));
  Energy kx = mz*theParent->sudPx() - cphi*pt;
  Energy ky = mz*theParent->sudPy() - sphi*pt; 
  double alpha = theParent->sudAlpha();

  // on-shell child
  c2->sudAlpha(mz*alpha);
  c2->sudPx(kx);
  c2->sudPy(ky);
  c2->sudBeta( (sqr(c2->data().mass())+c2->sudPperp2())/alpha/p_dot_n() );
  c2->set5Momentum(sudakov2Momentum(c2->sudAlpha(), c2->sudBeta(), 
				    c2->sudPx(), c2->sudPy()));
  // spacelike child:
  c1->sudAlpha(theParent->sudAlpha() - c2->sudAlpha());
  c1->sudBeta(theParent->sudBeta() - c2->sudBeta());
  c1->sudPx(theParent->sudPx() - c2->sudPx());
  c1->sudPx(theParent->sudPy() - c2->sudPy());
  c1->set5Momentum(theParent->momentum() - c2->momentum()); 

}


void IS_QtildaShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast ) {
  theLast->sudAlpha(1.0);
  theLast->sudBeta(0.0);
  theLast->sudPx(0.0*GeV);
  theLast->sudPy(0.0*GeV);
  theLast->set5Momentum(_pVector);
}


Energy IS_QtildaShowerKinematics1to2::jetMass() {

  Energy mass = Energy();

  //***LOOKHERE*** WRITE THE CODE 

  return mass;

}


