// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IS_QtildaShowerKinematics1to2 class.
//

#include "IS_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"

using namespace Herwig;

void IS_QtildaShowerKinematics1to2::
updateChildren( const tShowerParticlePtr theParent, 
		const ShowerParticleVector theChildren ) const {
  // set up initial properties
  if(theParent->children().empty()) {
    // no z for angular ordering in backward branchings
    const ShowerIndex::InteractionType inter=splittingFn()->interactionType();
    theParent->setEvolutionScale(inter, qtilde());
    theChildren[1]->setEvolutionScale(inter, (1.-z())*qtilde());
    // set proper colour connections
    splittingFn()->colourConnection(theParent,theChildren[0],theChildren[1],true);
    // set proper parent/child relationships
    theParent->addChild(theChildren[0]);
    theParent->addChild(theChildren[1]);
    theParent->x(theChildren[0]->x()/z());
    // Now fix the hadrons connections
    tPPtr hadron;
    if(theChildren[0]->parents().size() == 2) hadron = theChildren[0]->parents()[0];
    else throw Exception() << "IS_QtildaShowerKinematics1to2::"
			   << "updateChildren not one parent!" 
			   << Exception::runerror;
    hadron->abandonChild(theChildren[0]);
    hadron->addChild(theParent);
  }
  // update properties of children needed for branching
  // time-like child
  else {
    theChildren[1]->sudAlpha((1.-z())*theParent->sudAlpha());
    double cphi = cos(phi()),sphi = sin(phi());
    theChildren[1]->sudPx((1.-z())*theParent->sudPx() - cphi*pT());
    theChildren[1]->sudPy((1.-z())*theParent->sudPy() - sphi*pT());
    // space-like child
    theChildren[0]->sudAlpha(theParent->sudAlpha() - theChildren[1]->sudAlpha());
    theChildren[0]->sudBeta( theParent->sudBeta()  - theChildren[1]->sudBeta() );
    theChildren[0]->sudPx(   theParent->sudPx()    - theChildren[1]->sudPx()   );
    theChildren[0]->sudPy(   theParent->sudPy()    - theChildren[1]->sudPy()   );
  }
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
  double beta = 0.5*(sqr(c2->data().constituentMass())+c2->sudPperp2())
    /(c2->sudAlpha()*p_dot_n());
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
