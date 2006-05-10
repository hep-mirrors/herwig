// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "IS_QtildaShowerKinematics1to2.tcc"
#endif
#include "ShowerParticle.h"

using namespace Herwig;

IS_QtildaShowerKinematics1to2::~IS_QtildaShowerKinematics1to2() {}

void IS_QtildaShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector theChildren ) {
  // this is empty...
}


void IS_QtildaShowerKinematics1to2::
updateParent(const tShowerParticlePtr theParent, 
	     const ParticleVector theChildren ) {

  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
  // get shower variables from 1st child in order to keep notation
  // parent->(c1, c2) clean even though the splitting was initiated
  // from c1.  The name updateParent is still referring to the
  // timelike branching though.
  bool glueEmits = theParent->data().id() == ParticleID::g;
  Energy kinCutoff;
  if (glueEmits) {
    kinCutoff = showerVariables()->kinematicCutOff(kinScale(),c1->data().mass());
  } else {
    kinCutoff = showerVariables()->kinematicCutOff(kinScale(),theParent->data().mass());
  }
  double mz = 1.-c1->showerKinematics()->z();
  double cphi = cos(c1->showerKinematics()->phi());
  double sphi = sin(c1->showerKinematics()->phi());
  Energy pt = sqrt(sqr(mz*c1->showerKinematics()->qtilde()) - 
 		   (1.-mz)*sqr(kinCutoff));
  Energy kx = mz*theParent->sudPx() - cphi*pt;
  Energy ky = mz*theParent->sudPy() - sphi*pt; 
  double alpha = theParent->sudAlpha();
  // on-shell child
  c2->sudAlpha(mz*alpha);
  c2->sudPx(kx);
  c2->sudPy(ky);
  double beta = 0.0;
  if (c2->data().id() == ParticleID::g) {
    beta = (sqr(showerVariables()->gluonMass())+c2->sudPperp2())
      /c2->sudAlpha()/p_dot_n()/2.;
  } else {
    beta = (sqr(c2->data().constituentMass())+c2->sudPperp2())
      /c2->sudAlpha()/p_dot_n()/2.;
  }
  
  c2->sudBeta(beta);
  c2->set5Momentum(sudakov2Momentum(c2->sudAlpha(), c2->sudBeta(), 
				    c2->sudPx(), c2->sudPy()));
  // spacelike child:
  c1->sudAlpha(theParent->sudAlpha() - c2->sudAlpha());
  c1->sudBeta(theParent->sudBeta() - c2->sudBeta());
  c1->sudPx(theParent->sudPx() - c2->sudPx());
  c1->sudPy(theParent->sudPy() - c2->sudPy());
  Lorentz5Momentum pc1(theParent->momentum() - c2->momentum());
  pc1.rescaleMass();
  c1->set5Momentum(pc1);
}

void IS_QtildaShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast ) {
  if(theLast->isFinalState()) return;
  theLast->sudAlpha(theLast->x());
  theLast->sudBeta(sqr(theLast->data().mass())/
		   theLast->sudAlpha()/p_dot_n()/2.);
  theLast->sudPx(0.0*GeV);
  theLast->sudPy(0.0*GeV);
  theLast->set5Momentum(
    sudakov2Momentum(theLast->sudAlpha(), theLast->sudBeta(), 
		     theLast->sudPx(), theLast->sudPy()));
}

void IS_QtildaShowerKinematics1to2::calculatepT(vector<Energy>)
{
  throw Exception() << "IS_QtildaShowerKinematics1to2::calculatepT"
		    << "not yet implemented " << Exception::runerror;
}
