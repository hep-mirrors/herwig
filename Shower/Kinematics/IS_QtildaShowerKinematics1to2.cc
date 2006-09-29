// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "IS_QtildaShowerKinematics1to2.tcc"
#endif
#include "ShowerParticle.h"

using namespace Herwig;

void IS_QtildaShowerKinematics1to2::
updateChildren( const tShowerParticlePtr, 
		const ShowerParticleVector) const {
  // not needed
}


void IS_QtildaShowerKinematics1to2::
updateParent(const tShowerParticlePtr theParent, 
	     const ParticleVector theChildren ) const {

  ShowerParticlePtr c1 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[0]);
  ShowerParticlePtr c2 = dynamic_ptr_cast<ShowerParticlePtr>(theChildren[1]);
  // get shower variables from 1st child in order to keep notation
  // parent->(c1, c2) clean even though the splitting was initiated
  // from c1.  The name updateParent is still referring to the
  // timelike branching though.
  // on-shell child
  double beta = 0.0;
  if (c2->data().id() == ParticleID::g) 
    {beta = (sqr(showerVariables()->gluonMass())+c2->sudPperp2())
	/c2->sudAlpha()/p_dot_n()/2.;}
  else 
    {beta = (sqr(c2->data().constituentMass())+c2->sudPperp2())
      /c2->sudAlpha()/p_dot_n()/2.;}
  c2->sudBeta(beta);
  c2->set5Momentum(sudakov2Momentum(c2->sudAlpha(), c2->sudBeta(), 
				    c2->sudPx()   , c2->sudPy(),0));
  // spacelike child
  Lorentz5Momentum pc1(theParent->momentum() - c2->momentum());
  pc1.rescaleMass();
  c1->set5Momentum(pc1);
}

void IS_QtildaShowerKinematics1to2::
updateLast( const tShowerParticlePtr theLast,unsigned int iopt ) const {
  if(theLast->isFinalState()) return;
  theLast->sudAlpha(theLast->x());
  theLast->sudBeta(sqr(theLast->data().mass())/
		   theLast->sudAlpha()/p_dot_n()/2.);
  theLast->sudPx(0.0*GeV);
  theLast->sudPy(0.0*GeV);
  theLast->set5Momentum(sudakov2Momentum(theLast->sudAlpha(), theLast->sudBeta(), 
					 theLast->sudPx(), theLast->sudPy(),iopt));
}
