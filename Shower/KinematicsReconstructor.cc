// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KinematicsReconstructor class.
//

#include "KinematicsReconstructor.h"
#include "ShowerParticle.h"
#include "ShowerConfig.h"
#include "Pythia7/Interface/ClassDocumentation.h"
#include "Pythia7/Persistency/PersistentOStream.h"
#include "Pythia7/Persistency/PersistentIStream.h"
// #include "Pythia7/Interface/Parameter.h" 
#include "Pythia7/Interface/Reference.h" 
#include "Pythia7/Interface/RefVector.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "Pythia7/CLHEPWrap/Lorentz5Vector.h"

using namespace Herwig;


KinematicsReconstructor::~KinematicsReconstructor() {}


void KinematicsReconstructor::persistentOutput(PersistentOStream & os) const {
  os << _specialProcesses;  
}


void KinematicsReconstructor::persistentInput(PersistentIStream & is, int) {
  is >> _specialProcesses;
}


ClassDescription<KinematicsReconstructor> KinematicsReconstructor::initKinematicsReconstructor;
// Definition of the static class description member.


void KinematicsReconstructor::Init() {

  static ClassDocumentation<KinematicsReconstructor> documentation
    ( "This class is responsible for the kinematics reconstruction of the showering, ",
      "including the kinematics reshuffling necessary to compensate for the recoil of ",
      "the emissions." );
  
  static RefVector<KinematicsReconstructor,MEBase> interfaceVecSpecialProcesses
    ( "VecSpecialProcesses",
      "The collection (vector) of special processes.",
      &KinematicsReconstructor::_specialProcesses, 0, false, false, true, false );

}

// ------------------------------------------------------------------------

void KinematicsReconstructor::reconstructHardJets( const MapShower & mapShowerHardJets,
						   const Lorentz5Momentum & pBeamHadron1,
						   const Lorentz5Momentum & pBeamHadron2,
                                                   const tcMEPtr specialHardProcess ) 
  throw (Veto, Stop, Exception) {

  //***LOOKHERE***  First, treat the initial state, as follows.
  //                -------------------------------------------
  //                Loop over the map, and for each incoming hard particle 
  //                with flag true
  //                  call  reconstructSpaceLikeJet(...)  . 
  //                Then, if such method has being called at least once, then 
  //                If ( not specialHardSubprocess ) Then 
  //                    call  solveOverallCMframeBoost(...)
  //                Else
  //                    check which specialHardSubprocess we have.
  //                    In the case of D.I.S. 
  //                       call  solveSpecialDIS_CMframeBoost(...)
  //
  //                Second, do the overall CM boost, as follows.
  //                --------------------------------------------
  //                If ( not specialHardSubprocess ) Then  
  //                  Boost every outgoing jets (that is jets originated from 
  //                  an outgoing particle from the hard subprocess, but not
  //                  the forward jets coming from the initial state radiation)
  //                  from the Lab -> to the hard subprocess frame, using
  //                      pHard1Initial + pHard2Initial
  //                  as momentum of that frame. Then boost back, using this
  //                  time the new hard subprocess frame (obtained from the
  //                  method  solveOverallCMframeBoost(...) )
  //                      pHard1Final + pHard2Final .   
  //                Else
  //                    check which specialHardSubprocess we have.
  //                    In the case of D.I.S. , boost only the outgoing jet 
  //                    (that is the one originated from the outgoing parton from 
  //                     the hard subprocess lepton + parton -> lepton' + parton')
  //                    from the Lab -> to the hard subprocess frame, using
  //                        pLepton + pHardInitial
  //                    as momentum of that frame. Then boost back, using this
  //                    time the new hard subprocess frame (obtained from the
  //                    method  solveSpecialDIS_CMframeBoost(...) )
  //                        pLepton + pHardFinal .   
  // 
  //                Third, treat the final state, as follows.
  //                -----------------------------------------
  //                Loop again over the map, and for each outgoing hard particle 
  //                with flag true
  //                  call  reconstructTimeLikeJet(...)
  //                Then, if such method has being called at least once, then 
  //                  call  solveKfactor(...)
  //                and then  
  //                  call  solveBoost(...)  
  //                and then deep boost the particle with the boost returned 
  //                by the latter method.
  //               
  //                If any of the above methods fails, throw an Exception. 

}


void KinematicsReconstructor::reconstructDecayJets( const MapShower & mapShowerDecayJets )
  throw (Veto, Stop, Exception) {

  //***LOOKHERE***  Loop over the map, until you get the decaying particle. 
  //                Then if this has flag true
  //                  call  reconstructSpecialTimeLikeDecayingJet(...)  . 
  //                Loop again over the map, and for each decay product particle 
  //                with flag true
  //                  call  reconstructTimeLikeJet(...)
  //                Then, if any of the above methods has being called at least 
  //                once, then 
  //                  call  solveKfactor(...)
  //                and then  
  //                  call  solveBoost(...)  
  //                and then deep boost the decay product particles with the 
  //                boost returned by the latter method.
  //                If any of the above methods fails, throw an Exception. 

}


bool KinematicsReconstructor::reconstructTimeLikeJet( const tShoParPtr particleJetParent ) {
  bool isOK = true;
  
  // create a vector with pointers to the relevant ShowerKinematics objects.
  tCollecShoKinPtr theChildren; 
  theChildren.clear();

  //***LOOKHERE*** Following the tree of children and grand-children of the 
  //               particleJetParent, find the "reconstruction fixed points"
  //               (see the method that bears such name in the ShowerParticle class),
  //               that is either childless or decaying particles.
  //               From these then go back, following the parent pointers, 
  //               filling the empty bits of the ShowerKinematics objects
  //               by enforcing the energy-momentum conservation in each
  //               vertex. Stop when get back to the  particleJetParent
  //               at which point the mass of the jet should be known.
  
  if( !particleJetParent->isReconstructionFixedPoint() ) {
    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( CollecShoParPtr::const_iterator cit = particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit ) {
      // reconstrunct again for any child:
      if ( !reconstructTimeLikeJet( *cit ) ) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	  //generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: " 
	  //		     << "failed!" << endl; 
	}
	isOK = false; 
      }
      // apparently, fill theChildren 
      theChildren.push_back( (*cit)->showerKinematics() ); 
    }
  } else {
    // it is a reconstruction fixpoint, ie kinematical data has to be available 
    // check this
    // if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
      // generator()->log() << particleJetParent->showerKinematics()
    //  -> 'check some suitable variables';
    // }
  }

  // recursion has reached an endpoint, ie we can reconstruct the
  // kinematics from the children.
  particleJetParent->showerKinematics()->updateParent( theChildren );

  return isOK;
}


bool KinematicsReconstructor::reconstructSpaceLikeJet( const tShoParPtr particleJetParent ) {
  bool isOK = true;

  //***LOOKHERE*** The methods should be basically identical to the 
  //               previous one (because we define "parent" the particle
  //               closer to the hard subprocess, therefore the iteration
  //               through the parent/children pointers is exactly the same
  //               regardless whether we are moving forward or backward
  //               in time: the movement is the same w.r.t the hard subprocess!).
  //               The only difference can be eventually in the ShowerKinematics
  //               objects...

  return isOK;
}
 

bool KinematicsReconstructor::
reconstructSpecialTimeLikeDecayingJet( const tShoParPtr particleJetParent ) {
  bool isOK = true;

  //***LOOKHERE*** Differently from the two previous methods, in this
  //               case we have to reconstruct the kinematics from
  //               particleJetParent and following the children 
  //               until we end up with a "reconstruction fixed point"
  //               (see the method that bears such name in the ShowerParticle class),
  //               that is either childless or decaying particles.

  return isOK;
}
 

double KinematicsReconstructor::solveKfactor( const Lorentz5Momentum & cmMomentum, 
					      const VecMomentaPtr & jetsMomentaPtr ) {
  double k = 0.0;

  // ***LOOKHERE*** Go into the rest frame of cmMomentum, where M is the
  //                rest mass, and then solve numerically  k  from the 
  //                algebric equation:
  //                   SUM_i=1,n sqrt( mjet_i^2 + k * vec{p_i}^2 ) = M   

  return k;
}


Vector3 KinematicsReconstructor::
solveBoost( const double k, const Lorentz5Momentum & momentum ) {
  double bx = 0.0, by = 0.0, bz = 0.0;
  
  // ***LOOKHERE*** Find the boost that transforms:
  //                   vec{momentum} ---> k*vec{momentum}
  //                Returns (0,0,0) in the case no solution exists.

  return Vector3(bx,by,bz);
}


bool KinematicsReconstructor::
solveOverallCMframeBoost( const Lorentz5Momentum & pBeamHadron1,
			  const Lorentz5Momentum & pBeamHadron2,
			  const Lorentz5Momentum & pBeamParton1,
			  const Lorentz5Momentum & pBeamParton2,
			  const Lorentz5Momentum & pHard1Initial,
			  const Lorentz5Momentum & pHard2Initial,
			  const Lorentz5Momentum & pHard1Intermediate,
			  const Lorentz5Momentum & pHard2Intermediate,
			  Lorentz5Momentum & pHard1Final,
			  Lorentz5Momentum & pHard2Final ) {

  bool isOK = true;

  //***LOOKHERE*** WRITE THE CODE
  //               Notice that the method does not return the boost we are
  //               looking for, but instead pHard1Final and pHard2Final,
  //               from which we get the new center of mass frame of the
  //               hard subprocess: (pHard1Final + pHard2Final).
  //               Solving the boost from the initial to the final hard 
  //               subprocess frame:
  //               (pHard1Initial + pHard2Initial) -> (pHard1Final + pHard2Final)
  //               would be difficult, but anyway unnecessary because
  //               all we need to do is to boost all outgoing jets from
  //               the initial to the final hard subprocess frame, and
  //               this can be easily obtained as follows. First, we boost
  //               the outgoing jets from the Lab to the initial hard
  //               subprocess frame, and to that we need only
  //                 (pHard1Initial + pHard2Initial)
  //               Then, we assume that these momenta w.r.t. the hard
  //               subprocess frame remain unchanged, the only change
  //               being instead the motion of this reference frame w.r.t.
  //               the Lab. Therefore, we boost "back" these momenta 
  //               from the subprocess frame to the Lab, using 
  //                 (pHard1Final + pHard2Final)
 
  return isOK;
}


bool KinematicsReconstructor::
solveSpecialDIS_CMframeBoost( const Lorentz5Momentum & pLepton,
			      const Lorentz5Momentum & pBeamHadron,
			      const Lorentz5Momentum & pBeamParton,
			      const Lorentz5Momentum & pHardInitial,
			      const Lorentz5Momentum & pHardIntermediate,
			      Lorentz5Momentum & pHardFinal ) {
  
  bool isOK = true;

  //***LOOKHERE*** WRITE THE CODE
 
  return isOK;

}
