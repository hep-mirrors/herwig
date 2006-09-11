// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QtildaShowerKinematics1to2 class.
//

#include "Decay_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Base/ShowerParticle.h"


using namespace Herwig;

void 
Decay_QtildaShowerKinematics1to2::updateChildren(const tShowerParticlePtr theParent, 
					      const ShowerParticleVector theChildren ) const
{
  if(theChildren.size() != 2)
    throw Exception() <<  "Decay_QtildaShowerKinematics1to2::updateChildren() " 
 		      << "Warning! too many children!" << Exception::eventerror;
   // get the interaction type
  const ShowerIndex::InteractionType interaction =splittingFn()->interactionType();
  // copy scales etc
  Energy dqtilde = qtilde();
  double dz = z(); 
  double dphi = phi();
  theChildren[0]->setEvolutionScale(interaction, dqtilde);
  theChildren[1]->setEvolutionScale(interaction, (1.-dz)*dqtilde);
  theChildren[0]->setInitiatesTLS(false);
  theChildren[1]->setInitiatesTLS(false);
  // determine alphas of children according to interpretation of z
  theChildren[0]->sudAlpha( dz*theParent->sudAlpha() ); 
  theChildren[1]->sudAlpha( (1.-dz)*theParent->sudAlpha() ); 
  // set the values
  theChildren[0]->sudPx(   pT()*cos(dphi) +     dz *theParent->sudPx() );
  theChildren[0]->sudPy(   pT()*sin(dphi) +     dz *theParent->sudPy() );
  theChildren[1]->sudPx( - pT()*cos(dphi) + (1.-dz)*theParent->sudPx() );
  theChildren[1]->sudPy( - pT()*sin(dphi) + (1.-dz)*theParent->sudPy() );
}

void Decay_QtildaShowerKinematics1to2::
updateParent( const tShowerParticlePtr theParent, 
	      const ParticleVector theChildren ) const
{
  throw Exception() << "Decay_QtildaShowerKinematics1to2::updateParent not implemented"
		    << Exception::abortnow;
}

void Decay_QtildaShowerKinematics1to2::updateLast(const tShowerParticlePtr theLast,
						  unsigned int iopt) const
{
   // set beta component and consequently all missing data from that,
   // using the nominal (i.e. PDT) mass.
  Energy theMass = theLast->data().constituentMass(); 
  theLast->sudBeta( (sqr(theMass) + theLast->sudPperp2() 
   		     - sqr( theLast->sudAlpha() )*pVector().m2())
   		    / ( 2.*theLast->sudAlpha()*p_dot_n() ) );   
  // set that new momentum    
  theLast->set5Momentum(  sudakov2Momentum( theLast->sudAlpha(), theLast->sudBeta(), 
					    theLast->sudPx(), theLast->sudPy(),iopt));
}
