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
#include "Pythia7/Repository/EventGenerator.h"

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

bool KinematicsReconstructor::reconstructHardJets( const MapShower & mapShowerHardJets,
						   const Lorentz5Momentum & pBeamHadron1,
						   const Lorentz5Momentum & pBeamHadron2,
                                                   const tcMEPtr specialHardProcess ) 
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "KinematicsReconstructor::reconstructHardJets: "
		       << " ===> START DEBUGGING <=== " << endl;
  }
  
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
  // ***ACHTUNG!*** tried this below: 

  bool atLeastOnce = false;
  for ( MapShower::const_iterator cit = mapShowerHardJets.begin();
	cit != mapShowerHardJets.end(); ++cit ) {
    atLeastOnce = false;
    if ( cit->second && !cit->first->isFinalState() ) {
      atLeastOnce = ( atLeastOnce || reconstructSpaceLikeJet( cit->first ) ); 
    }
    if ( atLeastOnce ) {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	generator()->log() << "-- found initial state jets, try to boost them" << endl; 
      }      
      if ( !specialHardProcess ) {
	// solveOverallCMFrameBoost()
      } else { 
	// check which special hard process and call eg 
	// in the DIS case 
	// solveSpecialDIS_CMframeBoost()
      }
    }
  }

  //***LOOKHERE***  Second, do the overall CM boost, as follows.
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
  // ***ACHTUNG!*** hasn't been implemented since nothing special is
  // expected for the LEP case.

  //***LOOKHERE!*** Third, treat the final state, as follows.
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
  
  // collection of pointers to initial hard particle and jet momenta
  // for final boosts
  // CVecMomentaPtr jetsMomentaPtr; 
  //  VecMomenta initialMomenta; 
  Lorentz5Momentum p_cm = Lorentz5Momentum(); 
  JetKinVect jetKinematics;
  Lorentz5Momentum dum = Lorentz5Momentum();
  bool gottaBoost = false; 
  // only for debugging:
  Energy sum_qi = Energy(); 

  // find out whether we're in cm or not:
  for ( MapShower::const_iterator cit = mapShowerHardJets.begin();
	cit != mapShowerHardJets.end(); ++cit ) {
    p_cm += cit->first->momentum(); 
  }

  Vector3 beta_cm = p_cm.findBoostToCM();
  if ( beta_cm.mag() > 1e-12 ) gottaBoost = true;   

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "  p_cm = " << p_cm
		       << ", beta_cm = " << beta_cm
		       << ", boost? " << (gottaBoost ? "yes" : "no")
		       << endl; 
  }

  atLeastOnce = false;
  
  for ( MapShower::const_iterator cit = mapShowerHardJets.begin();
	cit != mapShowerHardJets.end(); ++cit ) {
    
    if ( cit->first->isFinalState() ) {            

      JetKinStruct tempJetKin;      
      tempJetKin.parent = cit->first; 
      tempJetKin.p = cit->first->momentum();
      if ( gottaBoost ) tempJetKin.p.boost( beta_cm ); 
      if ( cit->second ) 
	atLeastOnce = ( reconstructTimeLikeJet( cit->first ) || atLeastOnce ); 
      if ( gottaBoost ) {
	dum = cit->first->momentum();
	dum.boost( beta_cm ); 
	cit->first->momentum( dum );
      }
      tempJetKin.q = cit->first->momentum();       
      jetKinematics.push_back( tempJetKin );  

      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	sum_qi += tempJetKin.q.mag();
	generator()->log() << "  reconstructed "
			   << cit->first->data().PDGName()
			   << "-jet, q = "
			   << cit->first->momentum()
			   << endl
			   << "  remind momenta of "
			   << tempJetKin.parent
			   << ", should be " 
			   << cit->first 
			   << "."
			   << endl 
			   << "  p_i = " 
			   << tempJetKin.p
			   << ", q_i = " 
			   << tempJetKin.q
			   << endl; 
      }
    }
  }
 
  double k = 0.0; 
  if ( atLeastOnce ) {
    k = solveKfactor( p_cm.mag(), jetKinematics );
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  reshuffling with k = " << k << endl; 
      if ( k < 0. || k > 1. ) {
	generator()->log() << "  Warning! k outside 0..1. sum_qi > root_s ? " 
			   << (sum_qi > p_cm.mag() ? "yes" : "no") << endl
			   << "  VETO on this reconstruction -> restart shower! "
			   << endl;
      }
    }
    if ( HERWIG_DEBUG_LEVEL == HwDebug::minimal_Shower   
	 && generator()->currentEventNumber() < 1000 ) {    
      generator()->log() << "reshuffling with k = " << k; 
      if ( k < 0. || k > 1. ) {
	generator()->log() << " VETO -> restart shower! " << endl;
      }
    }
    if ( k < 0. || k > 1. ) return false; 
  }

  for ( JetKinVect::iterator it = jetKinematics.begin();
	it != jetKinematics.end(); ++it ) {
    LorentzRotation Trafo = LorentzRotation(); 
    if ( atLeastOnce ) Trafo = solveBoost(k, it->q, it->p);
    if ( gottaBoost ) Trafo.boost( -beta_cm );
    if ( atLeastOnce || gottaBoost ) it->parent->deepTransform( Trafo );
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    Lorentz5Momentum p_cm_after = Lorentz5Momentum(); 
    for ( MapShower::const_iterator cit = mapShowerHardJets.begin();
	  cit != mapShowerHardJets.end(); ++cit ) {
      p_cm_after += cit->first->momentum(); 
    }
    generator()->log() << "  reshuffling finished: p_cm  = " << p_cm << endl
		       << "                        p_cm' = " << p_cm_after << endl;
    p_cm_after = p_cm - p_cm_after; 
    if( sqr(p_cm_after.x()/MeV) > 1e-4 
	|| sqr(p_cm_after.y()/MeV) > 1e-4
	|| sqr(p_cm_after.z()/MeV) > 1e-4 
	|| sqr(p_cm_after.t()/MeV) > 1e-4 ) {
      generator()->log() << "  Warning! momentum conservation?!" 
			 << endl;
    } else {
      generator()->log() << "  ok!" 
			 << endl;      
    }
    generator()->log() << "KinematicsReconstructor::reconstructHardJets: "
		       << " ===> END DEBUGGING <=== " << endl;
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower 
       && generator()->currentEventNumber() < 1000) {    
    generator()->log() << ", p_cm = " << p_cm << endl;
    for ( JetKinVect::const_iterator it = jetKinematics.begin();
	  it != jetKinematics.end(); ++it ) {
      tCollecShoParPtr fs = it->parent->getFSChildren();
      generator()->log() << (it->parent)->data().PDGName()
			 << "-jet, Q = " << (it->parent)->momentum().m()
			 << ", q = " << (it->parent)->momentum()
 			 << ". " << fs.size() 
 			 << " ch:" << endl; 
      for ( tCollecShoParPtr::const_iterator jt = fs.begin(); 
	    jt != fs.end(); ++jt ) {
	generator()->log() << "  " << (*jt)->data().PDGName()
			   << ", q = " << (*jt)->momentum()
			   << endl; 
      }
    }
  }

  return true; 
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
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: "
		       << " ===> START DEBUGGING <=== " << endl;
  }

  //***LOOKHERE*** Following the tree of children and grand-children of the 
  //               particleJetParent, find the "reconstruction fixed points"
  //               (see the method that bears such name in the ShowerParticle class),
  //               that is either childless or decaying particles.
  //               From these then go back, following the parent pointers, 
  //               filling the empty bits of the ShowerKinematics objects
  //               by enforcing the energy-momentum conservation in each
  //               vertex. Stop when get back to the  particleJetParent
  //               at which point the mass of the jet should be known.
  
  if( !(particleJetParent->isReconstructionFixedPoint()) ) {

    // if not a reconstruction fixpoint, dig deeper for all children:
    for ( CollecShoParPtr::const_iterator cit = particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit ) {
      // reconstrunct again for any child:
      if ( !reconstructTimeLikeJet( *cit ) ) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	  generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: " 
	  		     << "failed!" << endl; 
	}
	isOK = false; 
      }
    }
  } else {

    // it is a reconstruction fixpoint, ie kinematical data has to be available 
    // check this    
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  * found fixpoint: "
			 << particleJetParent->data().PDGName()
			 << ", gets m = " 
			 << particleJetParent->data().mass()
			 << " #children = " 
			 << particleJetParent->children().size()
			 << endl;
    }

    // unfortunately the particleJetParent doesn't have his own
    // showerkinematics but this works fine and keeps the updateLast
    // in the ShowerKinematics class.
    particleJetParent->parent()->showerKinematics()->updateLast( particleJetParent );    
  }

  // recursion has reached an endpoint once, ie we can reconstruct the
  // kinematics from the children.
  if( !(particleJetParent->isReconstructionFixedPoint()) ) {

    particleJetParent->showerKinematics()
      ->updateParent( particleJetParent, particleJetParent->children() );      
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  * reconstruct "
			 << particleJetParent->data().PDGName()
			 << ", gets sqrt(q2) = " 
			 << particleJetParent->momentum().m()
			 << " #children = " 
			 << particleJetParent->children().size()
			 << endl;
    }
  }

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: "
			   << " ===> END DEBUGGING <=== " << endl;
  }

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
 


double KinematicsReconstructor::momConsEq(const double & k, const Energy & root_s, 
					  const JetKinVect & jets) {
  Energy dum = Energy(); 
  for( JetKinVect::const_iterator it = jets.begin(); it != jets.end(); ++it ) { 
    dum += sqrt( (it->q).mass2() + sqr(k)*((it->p).vect().mag2()) );
  }
  return( dum - root_s ); 
}


const double KinematicsReconstructor::solveKfactor( const Energy & root_s, 
						    const JetKinVect & jets) {
  Energy2 s = sqr(root_s);

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
    generator()->log() << "KinematicsReconstructor::solveKFactor: "
		       << "==> start debugging <== " << endl;
  }

  if ( jets.size() < 2) { 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	generator()->log() << "  Warning! called with < 2 jets!" 
			   << endl;
    }    
  } else if ( jets.size() == 2 ) {
    
    if ( momConsEq( 0.0, root_s, jets ) < 0.0 ) {
      if ( sqr((root_s - jets[0].p.t() - jets[1].p.t())/MeV) < 1.e-4
	   && sqr((jets[0].p.x()+jets[1].p.x())/MeV) < 1.e-4
	   && sqr((jets[0].p.y()+jets[1].p.y())/MeV) < 1.e-4
	   && sqr((jets[0].p.z()+jets[1].p.z())/MeV) < 1.e-4 ) {
	return sqrt( ( sqr(s - jets[0].q.m2() - jets[1].q.m2()) 
		       - 4.*jets[0].q.m2()*jets[1].q.m2() )
		     /(4.*s*jets[0].p.vect().mag2()) );      
      } else {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	  generator()->log() << "  Warning! 2 jets and not in cm-frame!" 
			     << endl
			     << "  root_s = " << root_s
			     << ", p0t = " << jets[0].p.t()
			     << ", p0t = " << jets[1].p.t()
			     << ", p0t+p1t = " << jets[0].p.t()+jets[1].p.t()
			     << endl
			     << "  (dE2, dpx2, dpy2, dpz2) = ("
			     << sqr((root_s - jets[0].p.t() - jets[1].p.t())/eV) << ", "
			     << sqr((jets[0].p.x()+jets[1].p.x())/eV) << ", "
			     << sqr((jets[0].p.y()+jets[1].p.y())/eV) << ", "
			     << sqr((jets[0].p.z()+jets[1].p.z())/eV) << ")" 
			     << endl;
	}
	return 0.0; 
      }
    } else {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	generator()->log() << "  Warning! can't find a k! return -1!" << endl 
			   << "KinematicsReconstructor::solveKFactor: "
			   << "==> end debugging <== " << endl;
      }
      return -1.; 
    }

  } else { // i.e. jets.size() > 2, numerically

    double k1 = 0.; 
    double k2 = 1.; 
    double k = 0.; 

    if ( momConsEq( k1, root_s, jets ) < 0.0 ) {
      while ( momConsEq( k2, root_s, jets ) < 0.0 ) {
	k1 = k2; 
	k2 *= 2;       
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	  generator()->log() << "  (k1, k2) = (" << k1 
			     << ", " << k2 << ") ...moved interval." 
			     << endl; 
	}
      }
      while ( fabs( (k1 - k2)/(k1 + k2) ) > 1.e-10 ) {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	  generator()->log() << "  (k1, k2) = (" << k1 
			     << ", " << k2 
			     << "), (eq1, eq2) = (" 
			     << momConsEq( k1, root_s, jets ) 
			     << ", "
			     << momConsEq( k2, root_s, jets ) 
			     << ")"
			     << endl; 
	} 
	if( momConsEq( k2, root_s, jets ) == 0. ) {
	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	    generator()->log() << "KinematicsReconstructor::solveKFactor: "
			       << "==> end debugging <== " << endl;
	  }
	  return k2; 
	} else {
	  k = (k1+k2)/2.;
	  if ( momConsEq( k, root_s, jets ) > 0 ) {
	    k2 = k;
	  } else {
	    k1 = k; 
	  } 
	}
      } 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	generator()->log() << "KinematicsReconstructor::solveKFactor: "
			   << "==> end debugging <== " << endl;
      }
      return k1; 	  
    } else {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	generator()->log() << "  Warning! haven't found a k! return -1!" << endl 
			   << "KinematicsReconstructor::solveKFactor: "
			   << "==> end debugging <== " << endl;
      }
      return -1.; 
    }
  }
  return -1.; 
}


Vector3 KinematicsReconstructor::
solveBoostBeta( const double k, const Lorentz5Momentum & newq, const Lorentz5Momentum & oldp ) {
  

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "KinematicsReconstructor::solveBoostBeta: "
		       << "==> start debugging <== " << endl;
  }

  // try something different, purely numerical first: 
  // a) boost to rest frame of newq, b) boost with kp/E

  Energy q = newq.vect().mag(); 
  Energy2 qs = sqr(q); 
  Energy2 Q2 = newq.m2(); 
  Energy kp = k*(oldp.vect().mag()); 
  Energy2 kps = sqr(kp); 

  // usually we take the minus sign, since this boost will be smaller.
  // we only require |k \vec p| = |\vec q'| which leaves the sign of
  // the boost open but the 'minus' solution gives a smaller boost
  // parameter, i.e. the result should be closest to the previous
  // result. this is to be changed if we would get many momentum
  // conservation violations at the end of the shower from a hard
  // process.
  double betam = (q*sqrt(qs + Q2) - kp*sqrt(kps + Q2))/(kps + qs + Q2); 

  // only for debug info:
  double betap = (q*sqrt(qs + Q2) + kp*sqrt(kps + Q2))/(kps + qs + Q2); 
  
  // move directly to 'return' 
  Vector3 beta = -betam*(k/kp)*oldp.vect();
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper. 

  // leave this out if it's running properly! 
  if ( betam >= 0 ) {

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  found beta (-, +) = ("
			 << betam << ", " << betap << ")" 
			 << endl
			 << "  directions (th, phi): p = ("
			 << oldp.theta() << ", " << oldp.phi() << ")" << endl
			 << "                        q = (" 
			 << newq.theta() << ", " << newq.phi() << ")" << endl;
    }
    
    Lorentz5Momentum test = newq;
    
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  - boosted     q' = "
			 << test.boost( beta ) << endl
			 << "  - constructed q' = ("
			 << k*oldp.vect().x()
			 << ","
			 << k*oldp.vect().y()
			 << ","
			 << k*oldp.vect().z()
			 << ";"
			 << sqrt(kps + Q2)
			 << ")" << endl
			 << "KinematicsReconstructor::solveBoostBeta: "
			 << "==> end debugging <== " << endl;
    }

    return beta;

  } else {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "  Warning! no beta found!"
			 << "KinematicsReconstructor::solveBoostBeta: "
			 << "==> end debugging <== " << endl;
    }

    return Vector3(0., 0., 0.); 
  }
}


LorentzRotation KinematicsReconstructor::
solveBoost( const double k, const Lorentz5Momentum & newq, const Lorentz5Momentum & oldp ) {
  
  Energy q = newq.vect().mag(); 
  Energy2 qs = sqr(q); 
  Energy2 Q2 = newq.m2(); 
  Energy kp = k*(oldp.vect().mag()); 
  Energy2 kps = sqr(kp); 

  double betam = (q*sqrt(qs + Q2) - kp*sqrt(kps + Q2))/(kps + qs + Q2); 
  Vector3 beta = -betam*(k/kp)*oldp.vect();
  // note that (k/kp)*oldp.vect() = oldp.vect()/oldp.vect().mag() but cheaper. 

  Hep3Vector ax = newq.vect().cross( oldp.vect() ); 
  double delta = newq.vect().angle( oldp.vect() );
  LorentzRotation R; 
  R.rotate( delta, ax ).boost( beta ); 

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    Lorentz5Momentum qprime( k*oldp.x(), k*oldp.y(), k*oldp.z(), 
			     sqrt(kps + Q2), sqrt(Q2) ); 
    Lorentz5Momentum test = newq;     
    qprime = qprime - (R*test);  

    generator()->log() << "KinematicsReconstructor::solveBoost: "
		       << "==> start debugging <== " << endl;
    if (k>1. || k<0) 
      generator()->log() << "  Warning! invalid k!"<< endl; 
    generator()->log() << "  Rotate around " << ax/ax.mag() 
		       << ", angle = " << delta << "." << endl
		       << "  constr-trans = " << qprime;
    if( sqr(qprime.x()/MeV) > 1e-4 
	|| sqr(qprime.y()/MeV) > 1e-4
	|| sqr(qprime.z()/MeV) > 1e-4 
	|| sqr(qprime.t()/MeV) > 1e-4 ) {
      generator()->log() << endl << "  Warning! constructed boost is inconsistent!" 
			 << endl;
    } else {
      generator()->log() << " ok!" 
			 << endl;      
    }    
    generator()->log() << "KinematicsReconstructor::solveBoost: "
		       << "==> end debugging <== " << endl;
  }
  return R;
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
