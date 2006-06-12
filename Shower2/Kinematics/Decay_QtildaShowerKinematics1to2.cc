// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Decay_QtildaShowerKinematics1to2 class.
//

#include "Decay_QtildaShowerKinematics1to2.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Shower2/SplittingFunctions/SplittingFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Decay_QtildaShowerKinematics1to2.tcc"
#endif
#include "ShowerParticle.h"


using namespace Herwig;

Decay_QtildaShowerKinematics1to2::~Decay_QtildaShowerKinematics1to2() {}



void 
Decay_QtildaShowerKinematics1to2::updateChildren(const tShowerParticlePtr theParent, 
					      const ShowerParticleVector theChildren ) 
{
  if(theChildren.size() != 2)
    throw Exception() <<  "Decay_QtildaShowerKinematics1to2::updateChildren() " 
 		      << "Warning! too many children!" << Exception::eventerror;
   // get the interaction type
  const ShowerIndex::InteractionType interaction = 
    theParent->splitFun()->interactionType();
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
  // determine transverse momenta of children
  Energy kinCutoff;
  vector<Energy> masses(3);
  kinCutoff = showerVariables()->kinematicCutOff(kinScale(),theParent->data().mass());
  // KMH - 15/05/06  - I think this should probably be checked...
  masses[0]=max(kinCutoff, theParent->data().mass());
  masses[1]=masses[0];
  masses[2]=kinCutoff;
  calculatepT(masses);
  // set the values
  theChildren[0]->sudPx(   pT()*cos(dphi) +     dz *theParent->sudPx() );
  theChildren[0]->sudPy(   pT()*sin(dphi) +     dz *theParent->sudPy() );
  theChildren[1]->sudPx( - pT()*cos(dphi) + (1.-dz)*theParent->sudPx() );
  theChildren[1]->sudPy( - pT()*sin(dphi) + (1.-dz)*theParent->sudPy() );
}

void Decay_QtildaShowerKinematics1to2::
updateParent( const tShowerParticlePtr theParent, 
	      const ParticleVector theChildren ) {
  throw Exception() << "Decay_QtildaShowerKinematics1to2::updateParent not implemented"
		    << Exception::abortnow;
}

void Decay_QtildaShowerKinematics1to2::updateLast(const tShowerParticlePtr theLast,
						  unsigned int iopt) {
   // set beta component and consequently all missing data from that,
   // using the nominal (i.e. PDT) mass.
  Energy theMass; 
  if ( theLast->data().id() == ParticleID::g ) theMass = showerVariables()->gluonMass();
  else theMass = theLast->data().constituentMass(); 
  theLast->sudBeta( (sqr(theMass) + theLast->sudPperp2() 
   		     - sqr( theLast->sudAlpha() )*pVector().m2())
   		    / ( 2.*theLast->sudAlpha()*p_dot_n() ) );   
  // set that new momentum    
  theLast->set5Momentum(  sudakov2Momentum( theLast->sudAlpha(), theLast->sudBeta(), 
					    theLast->sudPx(), theLast->sudPy(),iopt));
}

void Decay_QtildaShowerKinematics1to2::calculatepT(vector<Energy> masses)
{
  double dz=z(),domz(1.-dz);
  Energy2 pPerp2=sqr(qtilde()*domz)-sqr(masses[0]*domz)-sqr(masses[2])*dz;
  if(pPerp2<0.) throw Exception() << "Decay_QtildaShowerKinematics1to2::calculatepT()"
 				  << " Warning! Can't get p_perp, \n" 
 				  << "  z = " << dz << Exception::eventerror;
  pT(sqrt(pPerp2));
}
