// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KinematicsReconstructor class.
//

#include "KinematicsReconstructor.h"
#include "ShowerParticle.h"
#include "ShowerConfig.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
// #include "ThePEG/Interface/Parameter.h" 
#include "ThePEG/Interface/Reference.h" 
#include "ThePEG/Interface/RefVector.h"
#include "Herwig++/Utilities/HwDebug.h"
#include "ThePEG/CLHEPWrap/Lorentz5Vector.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

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


double myrap(const Lorentz5Momentum &p) {
  if (p.e() < p.z()) {
    return 999999999999.9;
  } else {
    return p.rapidity();
  }
}


void boostChain(tPPtr p, const Vector3 &bv) {
  if (p->parents()[0] && p->parents()[0]->id() == 21 
      || abs(p->parents()[0]->id()) < 7) {
    boostChain(p->parents()[0], bv);
  }
  p->boost(bv);
}


bool KinematicsReconstructor::
reconstructHardISJets(const MapShower &hardJets) 
  throw (Veto, Stop, Exception) {

  bool atLeastOnce = false;
  std::vector<Lorentz5Momentum> p;
  std::vector<Lorentz5Momentum> pq;
  std::vector<Lorentz5Momentum> p_in;
  MapShower::const_iterator cit;
  for(cit = hardJets.begin(); cit != hardJets.end(); cit++) {
    p_in.push_back(cit->first->momentum());
    atLeastOnce = false;
    cout << "reconstructHardJets..." << endl << flush;
    if(!cit->first->isFinalState()) {
      atLeastOnce |= reconstructSpaceLikeJet(cit->first);
      p.push_back(cit->first->momentum());
      if (cit->first->showerKinematics()) {
	pq.push_back(cit->first->showerKinematics()->getBasis()[0]);
      } else {
	//	pq.push_back(cit->first->momentum());
	if (cit->first->parents().size() > 0) {
	  pq.push_back(cit->first->parents()[0]->momentum());
	} else {
	  cerr << "  Shower/KinematicsReconstructor::reconstructHardJets: "
	       << "Warning, bad pq!!!"
	       << endl;
	  pq.push_back(cit->first->momentum());
	}
      }
    }
  }

  cout << "printing p..." << endl;
  for(unsigned int i = 0; i < p.size(); i++) cout << p[i] << endl;
  cout << "printing pq..." << endl;
  for(unsigned int i = 0; i < pq.size(); i++) cout << pq[i] << endl;
  cout << "printing p_in..." << endl;
  for(unsigned int i = 0; i < p_in.size(); i++) cout << p_in[i] << endl;

  cout << "  computing initial DY kinematics..." << endl;
  double x1, x2;
  x1 = p_in[0].z()/pq[0].z();
  x2 = p_in[1].z()/pq[1].z();

  Energy MDY = (p_in[0] + p_in[1]).m();
  Energy2 S = (pq[0]+pq[1]).m2();
  double yDY = (p_in[0] + p_in[1]).rapidity();

  cout << "  x1 = " << x1 << ", x2 = " << x2 << endl
       << "  MDY = " << MDY/GeV << " = " << sqrt(S*x1*x2)/GeV << endl
       << "  |yDY| = " << abs(yDY) << " = " << abs(0.5*log(x1/x2)) << endl;
  

  if(atLeastOnce && p.size() == 2 && pq.size() == 2) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "-- found initial state jets, try to boost them" 
			 << endl; 
    }      
    // find alphas and betas in terms of desired basis      
    Energy2 p12 = pq[0]*pq[1];
    //    Energy2 S = 2.*p12;
    double a1, a2, b1, b2;
    Lorentz5Momentum p1p, p2p;
    a1 = p[0]*pq[1]/p12;
    b1 = p[0]*pq[0]/p12;
    a2 = p[1]*pq[1]/p12;
    b2 = p[1]*pq[0]/p12;
    p1p = p[0] - a1*pq[0] - b1*pq[1];
    p2p = p[1] - a2*pq[0] - b2*pq[1];
    cout << "KinReco spacelike Jet setup reconstruction..." << endl
	 << "  p1  = " << p[0] << endl
	 << "      = " << a1 << " pq[0] + "
	 << b1 << " pq[1] " << endl 
	 << "        + " << p1p << endl 
	 << "      = " 
	 << a1*pq[0] + b1*pq[1] + p1p << endl
	 << "  p2  = " << p[1] << endl
	 << "      = " << a2 << " pq[0] + "
	 << b2 << " pq[1] " << endl 
	 << "        + " << p2p << endl 
	 << "      = " 
	 << a2*pq[0] + b2*pq[1] + p2p << endl
	 << "  Sudakov basis and decomposition checks..." << endl
	 << "  pq[0] = " << pq[0] << endl
	 << "  pq[1] = " << pq[1] << endl
	 << "  pq[0]^2 = " << pq[0]*pq[0] 
	 << ", pq[1]^2 = " << pq[1]*pq[1] 
	 << ", pq[0].pq[1] = " << p12 << endl
	 << "  p1p.pq[0] = " << p1p*pq[0]
	 << ", p1p.pq[1] = " << p1p*pq[1] << endl
	 << "  p2p.pq[0] = " << p2p*pq[0]
	 << ", p2p.pq[1] = " << p2p*pq[1] << endl
	 << "  p[0] + p[1] = " << p[0] + p[1] << endl;

    cout << "  We want to restore  M = " << MDY << endl;
    cout << "                      y = " << yDY << endl;
    cout << "  Now:                M = " << (p[0]+p[1]).m() << endl;
//     if ((p[0]+p[1]).m()/GeV < 0) {      
//       cout << "Throwing an exception..." << endl;      
//       throw Exception() << Exception::eventerror;
//     } else 
    cout << "                      y = " << myrap(p[0]+p[1]) << endl;
    // compute kappa12
    Energy2 A, B, C;
    double kp = 1.0, rad;
    A = a1*b2*S;
    B = sqr(MDY) - (a1*b1+a2*b2)*S - sqr(p1p+p2p);
    C = a2*b1*S; 
    rad = 1.-4.*A*C/sqr(B);
    if (rad >= 0) {
      kp = B/(2.*A)*(1.+sqrt(rad));
      cout << "  Resulting kp = " << kp << endl;
    } else {
      cout << "WARNING! Can't get kappa_pm!" << endl;
    }
    
    // now compute k1, k2
    double k1 = 1.0, k2 = 1.0;
    rad = kp*(b1+kp*b2)/(kp*a1+a2)*(x1/x2);   
    if (rad > 0) {
      k1 = sqrt(rad);
      k2 = kp/k1;
      cout << "  Plus:  k1 = " << k1 << ", k2 = " << k2 << endl; 
    } else {
      cout << "WARNING! Can't get k1p, k2p!" << endl;
    }

    double beta1 = getBeta((a1+b1), (a1-b1), 
			   (k1*a1+b1/k1), (k1*a1-b1/k1));
    double beta2 = getBeta((a2+b2), (a2-b2), 
			   (a2/k2+k2*b2), (a2/k2-k2*b2));
    cout << "  found boost parameters: beta1 = " << beta1 
	 << ", beta2 = " << beta2 << endl;
    
    if (pq[0].z() > 0) {beta1 = -beta1; beta2 = -beta2;}
    // check
    Lorentz5Momentum p0, p1;
    Vector3 betaboost;
    betaboost = Vector3(0, 0, beta1);
    p0 = p[0].boost(betaboost);
    betaboost = Vector3(0, 0, beta2);;
    p1 = p[1].boost(betaboost);
    cout << "  p0+p1 after boost..." << endl
	 << "  p0 = " << p0 << endl
	 << "  p1 = " << p1 << endl
	 << "  M/MDY = " << (p0+p1).m() 
	 << "/" << MDY 
	 << " = " << (p0+p1).m()/MDY << endl 
	 << "  y/yDY = " 
	 << myrap(p[0]+p[1])
	 << "/" << yDY 
      	 << " = " << myrap(p[0]+p[1])/yDY
	 << endl;
    cout << "  Check: ";
    if (abs((p0+p1).m()/MDY-1.0) > 1e-5) {
      cout << "M bad! delta = " << abs((p0+p1).m()/MDY-1.0) << ". ";
    } else {
      cout << "M ok. ";
    }
    if (abs(myrap(p0+p1)/yDY-1.0) > 1e-5) {
      cout << "y bad! delta = " << abs(myrap(p0+p1)/yDY-1.0) 
	   << "." << endl;
    } else {
      cout << "y ok." << endl;
    }

    tPVector toBoost;
    for(cit = hardJets.begin(); cit != hardJets.end(); cit++) {
      toBoost.push_back(cit->first);
    }

    // before boost
    cout << "Check before boosts.  toBoost contains last two partons. " << endl
	 << "  M = " 
	 << (toBoost[0]->momentum() + toBoost[1]->momentum()).m() 
	 << ", y = " 
	 << (toBoost[0]->momentum() + toBoost[1]->momentum()).rapidity()
	 << endl;

    cout << "toBoost[0]->id() = " 
	 << toBoost[0]->id()
	 << ", ->momentum() = "
	 << toBoost[0]->momentum() 
	 << endl
	 << "->children().size() = " 
	 << toBoost[0]->children().size()
	 << ", [0]->id() = " 
	 << toBoost[0]->children()[0]->id()
	 << ", [0]->momentum() = " 
	 << toBoost[0]->children()[0]->momentum()
	 << endl
	 << "toBoost[1]->id() = " 
	 << toBoost[1]->id()
	 << ", ->momentum() = " 
	 << toBoost[1]->momentum()
	 << endl
	 << "->children().size() = " 
	 << toBoost[1]->children().size()
	 << ", [0]->id() = " 
	 << toBoost[1]->children()[0]->id()
	 << ", momentum() = " 
	 << toBoost[1]->children()[0]->momentum()
	 << endl
	 << "total parton momentum = "
	 << toBoost[0]->momentum() + toBoost[1]->momentum() 
	 << endl;

    betaboost = Vector3(0, 0, beta1);
    boostChain(toBoost[0], betaboost);
    betaboost = Vector3(0, 0, beta2);
    boostChain(toBoost[1], betaboost);
    
    cout << "Check after boosts.  toBoost contains last two partons. " << endl
	 << "  M = " 
	 << (toBoost[0]->momentum() + toBoost[1]->momentum()).m() 
	 << ", y = " 
	 << (toBoost[0]->momentum() + toBoost[1]->momentum()).rapidity()
	 << endl;

    cout << "toBoost[0]->id() = " 
	 << toBoost[0]->id()
	 << ", ->momentum() = "
	 << toBoost[0]->momentum() 
	 << endl
	 << "->children().size() = " 
	 << toBoost[0]->children().size()
	 << ", [0]->id() = " 
	 << toBoost[0]->children()[0]->id()
	 << ", [0]->momentum() = " 
	 << toBoost[0]->children()[0]->momentum()
	 << endl
	 << "toBoost[1]->id() = " 
	 << toBoost[1]->id()
	 << ", ->momentum() = " 
	 << toBoost[1]->momentum()
	 << endl
	 << "->children().size() = " 
	 << toBoost[1]->children().size()
	 << ", [0]->id() = " 
	 << toBoost[1]->children()[0]->id()
	 << ", momentum() = " 
	 << toBoost[1]->children()[0]->momentum()
	 << endl
	 << "total parton momentum = " 
	 << toBoost[0]->momentum() + toBoost[1]->momentum() 
	 << endl;


    // consider DY pair...

    vector<Lorentz5Momentum> pDY;
    pDY.push_back(toBoost[0]->children()[0]->children()[0]->momentum());
    pDY.push_back(toBoost[0]->children()[0]->children()[1]->momentum());

    cout << "DY pair momenta, no boost: " << endl 
	 << "  p0 = " << pDY[0] << endl
	 << "  p1 = " << pDY[1] << endl
	 << "  p0+p1 = " << pDY[0] + pDY[1]  << endl
	 << "  M = " << (pDY[0] + pDY[1]).m() << endl
	 << "  y = " << (pDY[0] + pDY[1]).rapidity() << endl;
    
    Vector3 boostRest = (pDY[0] + pDY[1]).findBoostToCM();
    Vector3 boostNewF = (toBoost[0]->momentum() + toBoost[1]->momentum())
      .boostVector();
    (pDY[0].boost(boostRest)).boost(boostNewF);
    (pDY[1].boost(boostRest)).boost(boostNewF);
    
    cout << "DY pair momenta, after boost: " << endl 
	 << "  p0 = " << pDY[0] << endl
	 << "  p1 = " << pDY[1] << endl
	 << "  p0+p1 = " << pDY[0] + pDY[1]  << endl
	 << "  M = " << (pDY[0] + pDY[1]).m() << endl
	 << "  y = " << (pDY[0] + pDY[1]).rapidity() << endl;
    
    // actually boost DY Vector Boson and DY leptons:
    toBoost[0]->children()[0]->boost(boostRest);
    toBoost[0]->children()[0]->boost(boostNewF);
    toBoost[0]->children()[0]->children()[0]->boost(boostRest);
    toBoost[0]->children()[0]->children()[0]->boost(boostNewF);
    toBoost[0]->children()[0]->children()[1]->boost(boostRest);
    toBoost[0]->children()[0]->children()[1]->boost(boostNewF);    

    pDY.clear();
    pDY.push_back(toBoost[0]->children()[0]->children()[0]->momentum());
    pDY.push_back(toBoost[0]->children()[0]->children()[1]->momentum());

    cout << "DY pair momenta, after boost, check from Evt Record: " << endl 
	 << "  p0 = " << pDY[0] << endl
	 << "  p1 = " << pDY[1] << endl
	 << "  p0+p1 = " << pDY[0] + pDY[1]  << endl
	 << "  M = " << (pDY[0] + pDY[1]).m() << endl
	 << "  y = " << (pDY[0] + pDY[1]).rapidity() << endl;

    // 	 << pq[0] + pq[1] << endl
    // 	 << "               = pq[0] = " << pq0bb << endl
// 	 << "               + pq[1] = " << pq1bb << endl;
//     Lorentz5Momentum p1old, p1new, p2old, p2new;
//     Vector3 betaboost;
//     // for 'plus' solution
//     // first jet
//     p1old = a1*pq0bb + b1*pq1bb;
//     p1new = k1p*a1*pq0bb + b1/k1p*pq1bb; 
//     betaboost = -beta1p*p1old.vect()/p1old.vect().mag();
//         cout << "  betaboost = " << betaboost << endl
// 	 << "  p1old =      " << p1old  << endl
// 	 << "  p1new =      " << p1new << endl;
//     p1old.boost(betaboost);
//     cout << "  from boost = " << p1old << endl; 
//     p1new.boost(-betacm);
//     cout << "  p1 in original frame = " << p1new << endl;

//     // second jet
//     p2old = a2*pq0bb + b2*pq1bb;
//     p2new = k2p*a2*pq0bb + b2/k2p*pq1bb; 
//     betaboost = -beta2p*p2old.vect()/p2old.vect().mag();
//         cout << "  betaboost = " << betaboost << endl
// 	 << "  p2old =      " << p2old  << endl
// 	 << "  p2new =      " << p2new << endl;
//     p2old.boost(betaboost);
//     cout << "  from boost = " << p2old << endl;     
//     p2new.boost(-betacm);
//     cout << "  p2 in original frame = " << p2new << endl;

    
  }
  return true;
}


bool KinematicsReconstructor::reconstructHardJets(const MapShower &hardJets,
						  const Lorentz5Momentum &pB1,
						  const Lorentz5Momentum &pB2,
                                                  const tcMEPtr sHardProcess) 
  throw (Veto, Stop, Exception) {
  
  if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
    generator()->log() << "KinematicsReconstructor::reconstructHardJets: full _____________________________"<< endl; 
  }
  

  /********
   * Second, do the overall CM boost, as follows.
   * --------------------------------------------
   *  If not specialHardSubprocess Then  
   *    Boost every outgoing jet (that is jets originated from 
   *    an outgoing particle from the hard subprocess, but not
   *    the forward jets coming from the initial state radiation)
   *    from the Lab -> to the hard subprocess frame, using
   *        pHard1Initial + pHard2Initial
   *    as momentum of that frame. Then boost back, using this
   *    time the new hard subprocess frame (obtained from the
   *    method  solveOverallCMframeBoost(...) )
   *        pHard1Final + pHard2Final .   
   *  Else
   *      check which specialHardSubprocess we have.
   *      In the case of D.I.S., boost only the outgoing jet 
   *      (that is the one originated from the outgoing parton from 
   *       the hard subprocess lepton + parton -> lepton' + parton')
   *      from the Lab -> to the hard subprocess frame, using
   *          pLepton + pHardInitial
   *      as momentum of that frame. Then boost back, using this
   *      time the new hard subprocess frame (obtained from the
   *      method  solveSpecialDIS_CMframeBoost(...) )
   *          pLepton + pHardFinal .   
   *******/
  // ***ACHTUNG!*** hasn't been implemented since nothing special is
  // expected for the LEP case.

  /********
   * Third, treat the final state, as follows.
   * -----------------------------------------
   * Loop again over the map, and for each outgoing hard particle 
   * with flag true
   *   call  reconstructTimeLikeJet(...)
   * Then, if such method has being called at least once, then 
   *   call  solveKfactor(...)
   * and then  
   *   call  solveBoost(...)  
   * and then deep boost the particle with the boost returned 
   * by the latter method.
   * 
   * If any of the above methods fails, throw an Exception. 
   *******/ 

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
  MapShower::const_iterator cit;
  for(cit = hardJets.begin(); cit != hardJets.end(); ++cit) {
    p_cm += cit->first->momentum(); 
  }

  Vector3 beta_cm = p_cm.findBoostToCM();
  if(beta_cm.mag() > 1e-12) gottaBoost = true;   

  if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {
    generator()->log() << "  p_cm = " << p_cm
		       << ", beta_cm = " << beta_cm
		       << ", boost? " << (gottaBoost ? "yes" : "no")
		       << endl;
  }

  bool atLeastOnce = false;

  for(cit = hardJets.begin(); cit != hardJets.end(); cit++) {
    if(cit->first->isFinalState()) {
      JetKinStruct tempJetKin;      
      tempJetKin.parent = cit->first; 
      tempJetKin.p = cit->first->momentum();

      if(gottaBoost) tempJetKin.p.boost(beta_cm); 
      //atLeastOnce = (reconstructTimeLikeJet(cit->first) || atLeastOnce); 
      atLeastOnce |= reconstructTimeLikeJet(cit->first);
      tempJetKin.q = cit->first->momentum();       
      jetKinematics.push_back(tempJetKin);  

      if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {    
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
  if(atLeastOnce) {
    k = solveKfactor(p_cm.mag(), jetKinematics);
    if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {    
      generator()->log() << "  reshuffling with k = " << k << endl; 
      if(k < 0. || k > 1.) {
	generator()->log() << "  Warning! k outside 0..1. sum_qi > root_s ? " 
			   << (sum_qi > p_cm.mag() ? "yes" : "no") << endl
			   << "  VETO on this reconstruction -> restart shower! "
			   << endl;
      }
    }
    if(HERWIG_DEBUG_LEVEL == HwDebug::minimal_Shower &&  
       generator()->currentEventNumber() < 1000) {    
      generator()->log() << "reshuffling with k = " << k; 
      if(k < 0. || k > 1.) {
	generator()->log() << " VETO -> restart shower! " << endl;
      }
    }
    if(k < 0. || k > 1.) return false; 
  }

  for(JetKinVect::iterator it = jetKinematics.begin();
      it != jetKinematics.end(); ++it) {
    LorentzRotation Trafo = LorentzRotation(); 
    if(atLeastOnce) Trafo = solveBoost(k, it->q, it->p);
    if(gottaBoost) Trafo.boost(-beta_cm);
    if(atLeastOnce || gottaBoost) it->parent->deepTransform(Trafo);
  }

  if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {    
    Lorentz5Momentum p_cm_after = Lorentz5Momentum(); 
    for(cit = hardJets.begin(); cit != hardJets.end(); cit++) {
      p_cm_after += cit->first->momentum(); 
    }
    generator()->log() << "  reshuffling finished: p_cm  = " 
		       << p_cm << endl
		       << "                        p_cm' = " 
		       << p_cm_after << endl;
    p_cm_after = p_cm - p_cm_after; 
    if(sqr(p_cm_after.x()/MeV) > 1e-4 
	|| sqr(p_cm_after.y()/MeV) > 1e-4
	|| sqr(p_cm_after.z()/MeV) > 1e-4 
	|| sqr(p_cm_after.t()/MeV) > 1e-4 ) {
      generator()->log() << "  Warning! momentum conservation?!" 
			 << endl;
    } else {
      generator()->log() << "  ok!" 
			 << endl;      
    }
//     generator()->log() << "KinematicsReconstructor::reconstructHardJets: "
// 		       << " ===> END DEBUGGING <=== " << endl;
  }

  if(HERWIG_DEBUG_LEVEL >= HwDebug::minimal_Shower) {
    generator()->log() << ", p_cm = " << p_cm << endl;
    for ( JetKinVect::const_iterator it = jetKinematics.begin();
	  it != jetKinematics.end(); ++it ) {
      tShowerParticleVector fs = it->parent->getFSChildren();
      generator()->log() << (it->parent)->data().PDGName()
			 << "-jet, Q = " << (it->parent)->momentum().m()
			 << ", q = " << (it->parent)->momentum()
 			 << ". " << fs.size() 
 			 << " ch:" << endl; 
      Lorentz5Momentum sumch = Lorentz5Momentum(); 
      for ( tShowerParticleVector::const_iterator jt = fs.begin(); 
	    jt != fs.end(); ++jt ) {
	generator()->log() << "  " << (*jt)->data().PDGName()
			   << ", q = " << (*jt)->momentum()
			   << endl; 
	sumch += (*jt)->momentum(); 
      }
      generator()->log() << "  jet-parent q = " << sumch << endl
			 << "    children q = " 
			 << (it->parent)->momentum() << endl
			 << "          diff = " 
			 << (it->parent)->momentum() - sumch 
			 << endl;
    }
  }

  return true; 
}


void KinematicsReconstructor::reconstructDecayJets(const MapShower &decayJets)
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


bool KinematicsReconstructor::reconstructTimeLikeJet(const tShowerParticlePtr particleJetParent) {
  bool isOK = true;
  
//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
//     generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: "
// 		       << " ===> START DEBUGGING <=== " << endl;
//   }

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
    for ( ParticleVector::const_iterator cit = particleJetParent->children().begin();
	  cit != particleJetParent->children().end(); ++cit ) {
      // reconstrunct again for any child:
      if ( !reconstructTimeLikeJet( dynamic_ptr_cast<ShowerParticlePtr>(*cit)))
      {
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
      generator()->log() << "KinematicsReconstructor::recTLJet: ";
      generator()->log() << "fixpoint: "
			 << particleJetParent->data().PDGName()
			 << ", m = " 
			 << particleJetParent->momentum().mass()
			 << " #ch = " 
			 << particleJetParent->children().size()
			 << endl;
    }

    // unfortunately the particleJetParent doesn't have his own
    // showerkinematics but this works fine and keeps the updateLast
    // in the ShowerKinematics class.
    if (dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])) {      
//       cout << "b " 
// 	   << particleJetParent->id() << endl 
// 	   << particleJetParent->parents().size() << endl 
// 	   << "id = " << particleJetParent->parents()[0]->id() << endl 
// 	   << dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])
// 	   << endl 
// 	   << dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])->showerKinematics()      
// 	   << flush << endl; 
      if (dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])
	  ->showerKinematics()) {
	dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])
	  ->showerKinematics()->updateLast( particleJetParent );          
      }
    } else {
      Energy dm; 
      if (particleJetParent->id() == ParticleID::g) 
	// trying to get the effective gluon mass through with the 
	// usually unused 5th momentum component...
	dm = particleJetParent->momentum().mass();
      else
	dm = particleJetParent->data().constituentMass(); 
      if (dm != particleJetParent->momentum().m()) {
	Lorentz5Momentum dum =  particleJetParent->momentum();
	dum.setMass(dm); 
	dum.rescaleRho(); 
	particleJetParent->setMomentum(dum); 
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	  generator()->log() << "KinematicsReconstructor::recTLJet: ";
	  generator()->log() << "set untouched " 
			     << particleJetParent->data().PDGName()
			     << " to m = " 
			     << dm/GeV
			     << " GeV" << endl; 
	}
      } else {
	// *** ACHTUNG! *** find some way to tell the calling method
	// reconstructHardJets() that the particle doesn't need
	// rescaling, ie the subsequent 'rescaling/boost' could be
	// obsolete.  Only important when a hard process radiates
	// rarely.  Definitely not important for e+e- -> jets! 
	// simply returning 'false' collides with another check. 
      }
    }
  }

  // recursion has reached an endpoint once, ie we can reconstruct the
  // kinematics from the children.
  if( !(particleJetParent->isReconstructionFixedPoint()) ) {

    particleJetParent->showerKinematics()
      ->updateParent( particleJetParent, particleJetParent->children() );      
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "KinematicsReconstructor::recTLJet: ";
      generator()->log() << "recon "
			 << particleJetParent->data().PDGName()
			 << ", sqrt(q2) = " 
			 << particleJetParent->momentum().m()
			 << " #ch = " 
			 << particleJetParent->children().size()
			 << endl;
    }
  }

//   if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
// 	generator()->log() << "KinematicsReconstructor::reconstructTimeLikeJet: "
// 			   << " ===> END DEBUGGING <=== " << endl;
//   }

  return isOK;
}


bool KinematicsReconstructor::
reconstructSpaceLikeJet( const tShowerParticlePtr p) {
  bool isOK = true;
  if(abs(p->parents()[0]->id()) < 99) {
    cout << "recoSLJet part = " << p << ", going back..." << endl; 
    if (!reconstructSpaceLikeJet(dynamic_ptr_cast<ShowerParticlePtr>(
				 p->parents()[0]))) {
      isOK = false;
    }
    cout << "back done." << endl;
  } else {
    if(p->children().size() == 2) {
      cout << "recoSLJet part = " << p << ", last..." << endl;
      if (dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])) {
	dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])
	  ->showerKinematics()->updateLast(p);    
	cout << "last done." << endl << flush;
      } else {
	cout << "last done (no update, assuming particle hasn't split!)" 
	     << endl;
      }
    }
  }
  if(!p->isFromHardSubprocess() && 
     dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])) {
    cout << "  recoSLJet part = " << p << ", updating Children..." << endl;
    cout << "  " 
	 << p
	 << " (" << p->id() << ")"
	 << ", " 
	 << p->children()[0] 
	 << " (" << p->children()[0]->id() << ")"
	 << ", "
	 << p->children()[1] 
	 << " (" << p->children()[1]->id() << ")"
	 << endl;
    dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])
      ->showerKinematics()->updateParent(p, p->children());
    cout << "  p0    = " << p->momentum() << endl
	 << "  p1+p2 = " 
	 << p->children()[0]->momentum() + p->children()[1]->momentum()
	 << endl
	 << "  p1    = " 
	 << p->children()[0]->momentum() 
	 << endl 
	 << "  p2    = "
	 << p->children()[1]->momentum()
	 << endl;
    cout << "  update done." << endl;      
  }
  return isOK;
}
 

bool KinematicsReconstructor::
reconstructSpecialTimeLikeDecayingJet( const tShowerParticlePtr particleJetParent ) {
  bool isOK = true;

  //***LOOKHERE*** Differently from the two previous methods, in this
  //               case we have to reconstruct the kinematics from
  //               particleJetParent and following the children 
  //               until we end up with a "reconstruction fixed point"
  //               (see the method that bears such name in the ShowerParticle class),
  //               that is either childless or decaying particles.

  return isOK;
}
 


double KinematicsReconstructor::momConsEq(const double & k, 
					  const Energy & root_s, 
					  const JetKinVect & jets) {
  Energy dum = Energy(); 
  for(JetKinVect::const_iterator it = jets.begin(); it != jets.end(); ++it) 
    dum += sqrt( (it->q).m2() + sqr(k)*(it->p).vect().mag2() );
  return( dum - root_s ); 
}


const double KinematicsReconstructor::solveKfactor( const Energy & root_s, 
						    const JetKinVect & jets) {
  Energy2 s = sqr(root_s);

  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
    generator()->log() << "KinematicsReconstructor::solveKFactor: extreme ________________________________" << endl;
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
	generator()->log() << "  Warning! can't find a k! return -1!" << endl; 
// 			   << "KinematicsReconstructor::solveKFactor: "
// 			   << "==> end debugging <== " << endl;
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
// 	  if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
// 	    generator()->log() << "KinematicsReconstructor::solveKFactor: "
// 			       << "==> end debugging <== " << endl;
// 	  }
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
//       if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
// 	generator()->log() << "KinematicsReconstructor::solveKFactor: "
// 			   << "==> end debugging <== " << endl;
//       }
      return k1; 	  
    } else {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::extreme_Shower ) {    
	generator()->log() << "  Warning! haven't found a k! return -1!" << endl; 
// 			   << "KinematicsReconstructor::solveKFactor: "
// 			   << "==> end debugging <== " << endl;
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
    //    generator()->log() << "  boohoo!" << endl;
    generator()->log() << "KinematicsReconstructor::solveBoost full _______________________________________"<< endl; 
    if (k>1. || k<0) 
      generator()->log() << "  Warning! invalid k!"<< endl; 
    generator()->log() << "  Rotate around " << ((ax.mag()/MeV > 1e-4) ? ax/ax.mag() : ax)
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
//     generator()->log() << "KinematicsReconstructor::solveBoost: "
// 		       << "==> end debugging <== " << endl;
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
