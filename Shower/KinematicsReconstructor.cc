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
//#include <ofstream>
#include <cassert>

// for checks:
#include "Herwig++/Utilities/SmplHist.h"

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

namespace {
  double myrap(const Lorentz5Momentum &p) {
    if (p.e() < p.z()) {
      return 999999999999.9;
    } else {
      return p.rapidity();
    }
  }

// my new version, uses 'coloured()' but does the same as the old one.
  void boostChain(tPPtr p, const Vector3 &bv) {
    if (p->parents()[0] && p->parents()[0]->coloured()) 
      boostChain(p->parents()[0], bv);
    p->boost(bv);
  }


// Durham version
//   void boostChain(tPPtr p, const Vector3 &bv) {
//     if (p && p->coloured()) boostChain(p->parents()[0], bv);
//     if (p->children().size() == 2 && p->children()[0] && p->children()[1]) {
//       p->children()[0]->boost(bv);
//       p->children()[1]->boost(bv);
//     }
//   }

// initial version
// void boostChainOld(tPPtr p, const Vector3 &bv) {
//   if (p->parents()[0] && p->parents()[0]->id() == 21 
//       || abs(p->parents()[0]->id()) < 7) {
//     boostChain(p->parents()[0], bv);
//   }
//   p->boost(bv);
//}
}

bool KinematicsReconstructor::
reconstructHardISJets(const MapShower &hardJets) 
  throw (Veto, Stop, Exception) {
  // we need to catch everything else internally!
  // otherwise, we get std::unexpected() and std::abort()
  try {
    bool atLeastOnce = false;
    std::vector<Lorentz5Momentum> p, pq, p_in;
    MapShower::const_iterator cit;
    // loop over jet parents
    for(cit = hardJets.begin(); cit != hardJets.end(); ++cit) {
      p_in.push_back(cit->first->momentum());
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "KinematicsReconstructor::reconstructHardISJets..."
			   << endl;
      }
      if(!cit->first->isFinalState()) {
	atLeastOnce |= reconstructSpaceLikeJet(cit->first);
	p.push_back(cit->first->momentum());
	if (cit->first->showerKinematics()) {
	  pq.push_back(cit->first->showerKinematics()->getBasis()[0]);
	} else {
	  if (!cit->first->parents().empty()) {
	    //	    pq.push_back(cit->first->parents()[0]->momentum());
	    Energy etemp = cit->first->parents()[0]->momentum().pz();
	    Lorentz5Momentum ptemp = Lorentz5Momentum(0, 0, etemp, abs(etemp));
	    pq.push_back(ptemp);
	  } else {
	    generator()->log() 
	      << "Shower/KinematicsReconstructor::reconstructHardJets: "
	      << "Warning, bad pq!!!\n";
	    pq.push_back(cit->first->momentum());
	  }
	}
      }
    }

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  printing p..." << endl << flush;
      for(unsigned int i = 0; i < p.size(); i++) 
	generator()->log() << p[i] << endl;
      generator()->log() << "  printing pq..." << endl;
      for(unsigned int i = 0; i < pq.size(); i++) 
	generator()->log() << pq[i] << endl;
      generator()->log() << "  printing p_in..." << endl;
      for(unsigned int i = 0; i < p_in.size(); i++) 
	generator()->log() << p_in[i] << endl;
      generator()->log() << "  computing initial DY kinematics..." << endl;
    }

    double x1 = p_in[0].z()/pq[0].z();
    double x2 = p_in[1].z()/pq[1].z();

    Energy MDY = (p_in[0] + p_in[1]).m();
    Energy2 S = (pq[0]+pq[1]).m2();
    double yDY = myrap(p_in[0] + p_in[1]);

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  x1 = " << x1 << ", x2 = " << x2 << endl
			 << "  MDY = " << MDY/GeV << " = " << sqrt(S*x1*x2)/GeV 
			 << endl
			 << "  |yDY| = " << abs(yDY) << " = " 
			 << abs(0.5*log(x1/x2)) << endl;
    }

    if(! (atLeastOnce && p.size() == 2 && pq.size() == 2)) return true;


    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "-- found initial state jets, try to boost them" 
			 << endl; 
    }      
    // find alphas and betas in terms of desired basis      
    Energy2 p12 = pq[0]*pq[1];
    double a1 = p[0]*pq[1]/p12;
    double b1 = p[0]*pq[0]/p12;
    double a2 = p[1]*pq[1]/p12;
    double b2 = p[1]*pq[0]/p12;
    Lorentz5Momentum p1p = p[0] - a1*pq[0] - b1*pq[1];
    Lorentz5Momentum p2p = p[1] - a2*pq[0] - b2*pq[1];
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "p[1] = "
	//			 << setprecision(15)
			 << pq[1].x() << ", " 
			 << pq[1].y() << ", " 
			 << pq[1].z() << ", " 
			 << pq[1].t() << ", " 
			 << endl;
      generator()->log() << "KinReco spacelike Jet setup reconstruction..." 
			 << endl
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
			 << "  p[0] + p[1] = " << p[0] + p[1] << endl
			 << "  We want to restore  M = " << MDY << endl
			 << "                      y = " << yDY << endl
			 << "  Now:                M = " << (p[0]+p[1]).m() << endl;
    }
    if ((p[0]+p[1]).e()/GeV < 0) 
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() 
	  << "  E_DY/GeV = " << (p[0]+p[1]).e()/GeV 
	  << ", E_1 = " << p[0].e() 
	  << ", E_2 = " << p[1].e()
	  << " not returning false!" << endl;
      }
    //    if ((p[0]+p[1]).e()/GeV < 0) return false;
      
    //       cout << "Throwing an exception..." << endl;      
    //       throw Exception() << Exception::eventerror;
    //     } else 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      if (abs((p[0]+p[1]).pz()) > abs((p[0]+p[1]).e())) 
	generator()->log() << "  |pz| > |E|!!!" << endl;
      else
	generator()->log() 
	  << "                      y = " 
	  << myrap(p[0]+p[1]) << endl;
    }
    // compute kappa12
    // is this textbook method for solving a quadratic
    // numerically stable if 4AC ~= B^2 ? check Numerical Recipes
    double kp = 1.0;
    Energy2 A = a1*b2*S;
    Energy2 B = sqr(MDY) - (a1*b1+a2*b2)*S - sqr(p1p+p2p);
    Energy2 C = a2*b1*S; 
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "4AC = " << 4.*A*C << ", B^2 = " << B*B 
			 << ", 4AC - B^2 = " << 4.*A*C - B*B << endl;
    }
    double rad = 1.-4.*A*C/sqr(B);
    if (rad >= 0) {
      kp = B/(2.*A)*(1.+sqrt(rad));
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Resulting kp = " << kp << endl;
      }
    } else {
      cerr << "WARNING! Can't get kappa_pm!\n";
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "WARNING! Can't get kappa_pm!" << endl;
      }
    }
    
    // now compute k1, k2
    double k1 = 1.0, k2 = 1.0;
    rad = kp*(b1+kp*b2)/(kp*a1+a2)*(x1/x2);   
    if (rad > 0) {
      k1 = sqrt(rad);
      k2 = kp/k1;
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Plus:  k1 = " << k1 
			   << ", k2 = " << k2 << endl; 
      }
    } else {
      cerr << "  Plus:  k1 = " << k1 
	   << "WARNING! Can't get k1p, k2p!\n";
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
	generator()->log() << "  Plus:  k1 = " << k1 
			   << "WARNING! Can't get k1p, k2p!" << endl;
      }
    }

    double beta1 = getBeta((a1+b1), (a1-b1), 
			   (k1*a1+b1/k1), (k1*a1-b1/k1));
    double beta2 = getBeta((a2+b2), (a2-b2), 
			   (a2/k2+k2*b2), (a2/k2-k2*b2));
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  Plus:  k1 = " << k1 
			 << "  found boost parameters: beta1 = " << beta1 
			 << ", beta2 = " << beta2 << endl;
    }
    
    if (pq[0].z() > 0) {beta1 = -beta1; beta2 = -beta2;}
    // check
    Vector3 betaboost(0, 0, beta1);
    Lorentz5Momentum p0 = p[0].boost(betaboost);
    betaboost = Vector3(0, 0, beta2);
    Lorentz5Momentum p1 = p[1].boost(betaboost);
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  initial check: p0+p1 after boost..." << endl
			 << "  p0 = " << p0 << endl
			 << "  p1 = " << p1 << endl
			 << "  M/MDY = " << (p0+p1).m() 
			 << "/" << MDY 
			 << " = " << (p0+p1).m()/MDY << endl 
			 << "  y/yDY = " 
			 << myrap(p[0]+p[1])
			 << "/" << yDY 
			 << " = " << myrap(p[0]+p[1])/yDY
			 << endl
			 << "  Check: ";
      if (abs((p0+p1).m()/MDY-1.0) > 1e-5) {
	generator()->log() << "M bad! delta = " 
			   << abs((p0+p1).m()/MDY-1.0) << ". ";
      } else {
	generator()->log() << "M ok. ";
      }
      if (abs(myrap(p0+p1)/yDY-1.0) > 1e-5) {
	generator()->log() << "y bad! delta = " << abs(myrap(p0+p1)/yDY-1.0) 
			   << "." << endl << flush;
      } else {
	generator()->log() << "y ok." << endl;
      }
    }

    tPVector toBoost;
    for(cit = hardJets.begin(); cit != hardJets.end(); ++cit) {
      toBoost.push_back(cit->first);
    }

    // before boost
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() 
	<< "Check before boosts.  toBoost contains last two partons. " << endl
	<< "  M = " 
	<< (toBoost[0]->momentum() + toBoost[1]->momentum()).m();
      if (abs((toBoost[0]->momentum() + toBoost[1]->momentum()).e()) 
	  > abs((toBoost[0]->momentum() + toBoost[1]->momentum()).pz())) 
	generator()->log() 
	  << ", y = " 
	  << myrap(toBoost[0]->momentum() + toBoost[1]->momentum())
	  << endl;
      else       
	generator()->log() << "y non-calculable!" << endl;
    }
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "toBoost[0]->id() = " 
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
    }
    betaboost = Vector3(0, 0, beta1);
    boostChain(toBoost[0], betaboost);
    betaboost = Vector3(0, 0, beta2);
    boostChain(toBoost[1], betaboost);
    
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "Check after boosts.  toBoost contains last two partons. " << endl
			 << "  M = " 
			 << (toBoost[0]->momentum() + toBoost[1]->momentum()).m() 
			 << ", y = " 
			 << myrap(toBoost[0]->momentum() + toBoost[1]->momentum())
			 << endl;
    }
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "  toBoost[0]->id() = " 
			 << toBoost[0]->id()
			 << ", ->momentum() = "
			 << toBoost[0]->momentum() 
			 << endl
			 << "  ->children().size() = " 
			 << toBoost[0]->children().size()
			 << ", [0]->id() = " 
			 << toBoost[0]->children()[0]->id()
			 << endl 
			 << ", [0]->momentum() = " 
			 << toBoost[0]->children()[0]->momentum()
			 << endl
			 << "  toBoost[1]->id() = " 
			 << toBoost[1]->id()
			 << ", ->momentum() = " 
			 << toBoost[1]->momentum()
			 << endl
			 << "  ->children().size() = " 
			 << toBoost[1]->children().size()
			 << ", [0]->id() = " 
			 << toBoost[1]->children()[0]->id()
			 << endl 
			 << ", momentum() = " 
			 << toBoost[1]->children()[0]->momentum()
			 << endl
			 << "total parton momentum = " 
			 << toBoost[0]->momentum() + toBoost[1]->momentum() 
			 << endl;
    }

    // consider DY pair...
    // this is specific to the default hard process
    vector<Lorentz5Momentum> pDY;
    pDY.push_back(toBoost[0]->children()[0]->children()[0]->momentum());
    pDY.push_back(toBoost[0]->children()[0]->children()[1]->momentum());

    // need to apply boosts to the colour singlet system
    Lorentz5Momentum psinglet(toBoost[0]->children()[0]->momentum());

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "DY pair momenta, no boost: " << endl 
			 << "  p0 = " << pDY[0] << endl
			 << "  p1 = " << pDY[1] << endl
			 << "  p0+p1 = " << pDY[0] + pDY[1]  << endl
			 << "  psinglet = " << psinglet << '\n'
			 << "  M = " << (pDY[0] + pDY[1]).m() << endl
			 << "  y = " << myrap(pDY[0] + pDY[1]) << endl;
    }
    Vector3 boostRest = psinglet.findBoostToCM();
    Vector3 boostNewF = (toBoost[0]->momentum() + toBoost[1]->momentum())
      .boostVector();
    
    // this is specific to the current process, see above
    (pDY[0].boost(boostRest)).boost(boostNewF);
    (pDY[1].boost(boostRest)).boost(boostNewF);
    
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "DY pair momenta, after boost: " << endl 
    			 << "  p0 = " << pDY[0] << endl
    			 << "  p1 = " << pDY[1] << endl
    			 << "  p0+p1 = " << pDY[0] + pDY[1]  << endl
    			 << "  M = " << (pDY[0] + pDY[1]).m() << endl
    			 << "  y = " << myrap(pDY[0] + pDY[1]) << endl;
    }
    
    // actually boost DY Vector Boson and DY leptons:
    toBoost[0]->children()[0]->boost(boostRest);
    toBoost[0]->children()[0]->boost(boostNewF);
    toBoost[0]->children()[0]->children()[0]->boost(boostRest);
    toBoost[0]->children()[0]->children()[0]->boost(boostNewF);
    toBoost[0]->children()[0]->children()[1]->boost(boostRest);
    toBoost[0]->children()[0]->children()[1]->boost(boostNewF);    
    // Changed by Durham group
    //    toBoost[0]->children()[0]->deepBoost(boostRest);
    //    toBoost[0]->children()[0]->deepBoost(boostNewF);
// #define PHILSCODE
#ifdef PHILSCODE
    // find remnants and repair their kinematics
    // PJS: Ignore the remnant for now, it must be built after the ISR

    // comment block one was here (see below)

    pDY.clear();
    pDY.push_back(toBoost[0]->children()[0]->children()[0]->momentum());
    pDY.push_back(toBoost[0]->children()[0]->children()[1]->momentum());

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() 
	<< "DY pair momenta, after boost, check from Evt Record: " << endl 
	<< "  p0 = " << pDY[0] << endl
	<< "  p1 = " << pDY[1] << endl
	<< "  p0+p1 = " << pDY[0] + pDY[1]  << endl
	<< "  M = " << (pDY[0] + pDY[1]).m() << endl
	<< "  y = " << myrap(pDY[0] + pDY[1]) << endl;
    }
    // repair remnant's momenta
    //     cout << "repairing 1st remnant..." << endl;
    //     toBoost[0]->parents()[0]->children()[0]
    //       ->setMomentum(toBoost[0]->parents()[0]->momentum()
    // 		    -toBoost[0]->momentum());
    //     cout << "repairing 2nd remnant..." << endl;
    //     toBoost[1]->parents()[0]->children()[0]
    //       ->setMomentum(toBoost[1]->parents()[0]->momentum()
    // 		    -toBoost[1]->momentum());

    ///////////////////////////////////////////////////////////////////////////////
#endif
#undef PHILSCODE

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() 
	<< "For checks we want to consider partons as well..." << endl
	<< toBoost[0]->id() << ", " << toBoost[0]->momentum() << endl
	<< toBoost[1]->id() << ", " << toBoost[1]->momentum() << endl
	<< "  sum = " << toBoost[0]->momentum()+toBoost[1]->momentum() 
	<< endl;
    } 
    
    tPPtr par1 = toBoost[0]->parents()[0];
    tPPtr par2 = toBoost[1]->parents()[0];
    int ia, ib; 
    if (par1->id() > 0) {ia = 0; ib = 1;}
    else {ia = 1; ib = 0;}

    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << "X, " << par1 << "(" << par1->id() << "), " 
			 << par2 << "(" << par2->id() << ")" << endl;
    }
    while(abs(par1->id()) < 99 && !par1->parents().empty()) 
      par1 = par1->parents()[0];
    while(abs(par2->id()) < 99 && !par2->parents().empty()) 
      par2 = par2->parents()[0];
    //if(!par1 || !par2) return true;
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      generator()->log() << par1->id() << ", " << par2->id() << endl;
      generator()->log() << par1->children()[0]->id() << ", ";
      if(par1->children().size() >= 2) 
	generator()->log() << par1->children()[1]->id() << ", ";
      else return true;
      generator()->log() << par2->children()[0]->id();
      if(par2->children().size() >= 2)
	generator()->log() << ", " << par2->children()[1]->id();
      else return true;
      generator()->log() << endl;
    }

    // comment block 2 was here

    Lorentz5Momentum q;
    double xp, xm;


    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
      // print whole emission chain
      generator()->log() << "Printing emission chains xp, xm, qperp, virt, 4mom (all in GeV)..." << endl;
      
      Lorentz5Momentum delta;
    
      generator()->log() << "leg 1:" << endl
			 << setprecision(4) 
			 <<   "H    " << par1->momentum()/GeV << endl; 
      q = par1->children()[0]->momentum();
      xp = q*pq[0]/p12; xm = q*pq[1]/p12;
      generator()->log() << "|=== " 
			 << setw(11) << xp 
			 << setw(11) << xm 
			 << setw(6) << q.perp()/GeV
			 << setw(7) << q.m()/GeV 
			 << setw(5) << par1->children()[0]->id()  
			 << " " << q/GeV << endl;
      q = par1->children()[1]->momentum();
      xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
      generator()->log() << "|    " 
			 << setw(11) << xp 
			 << setw(11) << xm 
			 << setw(6) << q.perp()/GeV
			 << setw(7) << q.m()/GeV 
			 << setw(5) << par1->children()[1]->id()  
			 << " " << q/GeV << endl;
      //     delta = par1->momentum() - par1->children()[0]->momentum() - par1->children()[1]->momentum();
      //     generator()->log() << "delta = " << delta/GeV << endl;
      par1 = par1->children()[1];
      //      generator()->log() << par1->children().size() <<endl;
      while (par1->children().size() == 2) {
	q = par1->children()[1]->momentum();
	xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
	generator()->log() << "|--- " 
			   << setw(11) << xp 
			   << setw(11) << xm 
			   << setw(6) << q.perp()/GeV
			   << setw(7) << q.m()/GeV 
			   << setw(5) << par1->children()[1]->id()  
			   << " " << q/GeV << endl;
	q = par1->children()[0]->momentum();
	xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
	generator()->log() << "|    " 
			   << setw(11) << xp 
			   << setw(11) << xm 
			   << setw(6) << q.perp()/GeV
			   << setw(7) << q.m()/GeV 
			   << setw(5) << par1->children()[0]->id()  
			   << " " << q/GeV << endl;
	//       delta = par1->momentum() - par1->children()[0]->momentum() - par1->children()[1]->momentum();
	//       generator()->log() << "delta = " << delta/GeV << endl;
	par1 = par1->children()[0];
      }
      
      generator()->log() << "leg 2:" << endl
			 <<   "H    " << par2->momentum()/GeV << endl;
      q = par2->children()[0]->momentum();
      xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
      generator()->log() << "|=== " 
			 << setw(11) << xp 
			 << setw(11) << xm 
			 << setw(6) << q.perp()/GeV
			 << setw(7) << q.m()/GeV 
			 << setw(5) << par2->children()[0]->id() 
			 << " " << q/GeV << endl;
      q = par2->children()[1]->momentum();
      xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
      generator()->log() << "|    " 
			 << setw(11) << xp 
			 << setw(11) << xm 
			 << setw(6) << q.perp()/GeV
			 << setw(7) << q.m()/GeV 
			 << setw(5) << par2->children()[1]->id()  
			 << " " << q/GeV << endl;
      //     delta = par2->momentum() - par2->children()[0]->momentum() 
      //       - par2->children()[1]->momentum();
      //     generator()->log() << "delta = " << delta/GeV << endl;
      par2 = par2->children()[1];
      while (par2->children().size() == 2) {
	q = par2->children()[1]->momentum();
	xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
	generator()->log() << "|--- " 
			   << setw(11) << xp 
			   << setw(11) << xm 
			   << setw(6) << q.perp()/GeV
			   << setw(7) << q.m()/GeV 
			   << setw(5) << par2->children()[1]->id()  
			   << " " << q/GeV << endl;
	q = par2->children()[0]->momentum();
	xp = q*pq[0]/p12; xm = q*pq[1]/p12; 
	generator()->log() << "|    " 
			   << setw(11) << xp 
			   << setw(11) << xm 
			   << setw(6) << q.perp()/GeV
			   << setw(7) << q.m()/GeV 
			   << setw(5) << par2->children()[0]->id()  
			   << " " << q/GeV << endl;
	//      delta = par2->momentum() - par2->children()[0]->momentum() 
	//	- par2->children()[1]->momentum();
	//      generator()->log() << "delta = " << delta/GeV << endl;
	par2 = par2->children()[0];
      }
    }

    // comment block 3 was here
  }
  catch(std::exception & e) {
    throw Exception() << "Caught exception\n"
		      << e.what() 
		      <<  "\nin KinematicsReconstructor::reconstructHardISJets"
		      << Exception::eventerror;
  }

  return true;
}



bool KinematicsReconstructor::reconstructHardJets(const MapShower &hardJets,
						  const Lorentz5Momentum &pB1,
						  const Lorentz5Momentum &pB2)
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

  // only for debugging:
  Energy sum_qi = Energy(); 

  // find out whether we're in cm or not:
  MapShower::const_iterator cit;
  Lorentz5Momentum p_cm = Lorentz5Momentum(); 
  for(cit = hardJets.begin(); cit != hardJets.end(); ++cit) {
    if (cit->first->isFinalState()) // avoids double counting if ISR and FSR are on
      p_cm += cit->first->momentum(); 
  }

  Vector3 beta_cm = p_cm.findBoostToCM();
  bool gottaBoost = (beta_cm.mag() > 1e-12);

  if(HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower) {
    generator()->log() << "  p_cm = " << p_cm
		       << ", beta_cm = " << beta_cm
		       << ", boost? " << (gottaBoost ? "yes" : "no")
		       << endl;
  }

  bool atLeastOnce = false;
  // collection of pointers to initial hard particle and jet momenta
  // for final boosts
  JetKinVect jetKinematics; 
  for(cit = hardJets.begin(); cit != hardJets.end(); cit++) {
    if(cit->first->isFinalState()) {
      JetKinStruct tempJetKin;      
      tempJetKin.parent = cit->first; 
      tempJetKin.p = cit->first->momentum();
      if(gottaBoost) tempJetKin.p.boost(beta_cm); 
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

  //  return true;
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
    throw Exception() << "Need to implemente KinematicsReconstructor"
		      << "::reconstructDecayJets()" << Exception::runerror;

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
    ShowerParticlePtr jetGrandParent = dynamic_ptr_cast<ShowerParticlePtr>
      (particleJetParent->parents()[0]);
    if (jetGrandParent) {      
//       cout << "b " 
// 	   << particleJetParent->id() << endl 
// 	   << particleJetParent->parents().size() << endl 
// 	   << "id = " << particleJetParent->parents()[0]->id() << endl 
// 	   << dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])
// 	   << endl 
// 	   << dynamic_ptr_cast<ShowerParticlePtr>(particleJetParent->parents()[0])->showerKinematics()      
// 	   << flush << endl; 
      if (jetGrandParent->showerKinematics()) {
	jetGrandParent->showerKinematics()->updateLast( particleJetParent );
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
  tShowerParticlePtr child;
  if(abs(p->parents()[0]->id()) < 99) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "recoSLJet part = " 
			 << p << ", going back to..." 
			 << *p->parents()[0] << endl; 
    }
    // NOTE: PJS - Added this check for final state showering
    tShowerParticlePtr parent = dynamic_ptr_cast<ShowerParticlePtr>
      (p->parents()[0]);
    if(parent) isOK = reconstructSpaceLikeJet(parent);
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log() << "back done." << endl;
    }
  } else {
    if(p->children().size() == 2) {
      if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	generator()->log() << "recoSLJet part = " << p 
			   << ", last..." << endl;
      }
      child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0]);
      if (child) {
	child->showerKinematics()->updateLast(p);    
// 	cerr << p->children()[0] << endl
// 	     << dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0])
// 	  ->showerKinematics() << endl;
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	  generator()->log() << "last done." << endl << flush;
	}
      } else {
	if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
	  generator()->log() 
	    << "last done (no update, assuming particle hasn't split!)" 
	    << endl;
	}
      }
    }
  }
  child = dynamic_ptr_cast<ShowerParticlePtr>(p->children()[0]);
  if(!p->isFromHardSubprocess() && child) {
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log()
	<< "  recoSLJet part = " 
	<< p << ", updating Children..." << endl
	<< "  " 
	<< p
	<< " (" << p->id() << ")"
	<< ", " 
	<< p->children()[0] 
	<< " (" << p->children()[0]->id() << ")"
	<< ", "
	<< p->children()[1] 
	<< " (" << p->children()[1]->id() << ")"
	<< endl;
    }
    child->showerKinematics()->updateParent(p, p->children());
    if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {    
      generator()->log()
	<< "  p0    = " << p->momentum() << endl
	<< "  p1+p2 = " 
	<< p->children()[0]->momentum() + p->children()[1]->momentum()
	<< endl
	<< "  p1    = " 
	<< p->children()[0]->momentum() 
	<< endl 
	<< "  p2    = "
	<< p->children()[1]->momentum()
	<< endl
	<< "  update done." << endl;      
    }
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

  throw Exception() << "need to implement "
		    << "KinematicsReconstructor::reconstructSpecialTimeLike"
		    << "DecayingJet()"
		    << Exception::runerror;
  

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
    return -1.0;
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
    // check convergence, if it's a problem maybe use Newton iteration?
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
  assert(k > 0.0 && k <= 1.0);
  
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
      
  throw Exception() << "need to implement "
		    << "KinematicsReconstructor::solveOverallCMframeBoost()"
		    << Exception::runerror;
  
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
  throw Exception() << "need to implement "
		    << "KinematicsReconstructor::solveSpecialDIS"
		    << "_CMframeBoost()"
		    << Exception::runerror;
  
  return isOK;

}

// comment block 1


    //     tPPtr r1, r2, q1, q2;
    //     r1 = toBoost[0]->parents()[0];
    //     r2 = toBoost[1]->parents()[0];
    //     while(abs(r1->id()) < 99) r1 = r1->parents()[0];
    //     while(abs(r2->id()) < 99) r2 = r2->parents()[0];
    //     if (r1->children()[0]->id() > 99) {
    //       r1 = r1->children()[0];
    //       q1 = r1->parents()[0]->children()[1];
    //     } else {
    //       r1 = r1->children()[1];
    //       q1 = r1->parents()[0]->children()[0];
    //     }
    //     if (r2->children()[0]->id() > 99) {
    //       r2 = r2->children()[0];
    //       q2 = r2->parents()[0]->children()[1];
    //     } else {
    //       r2 = r2->children()[1];
    //       q2 = r2->parents()[0]->children()[0];
    //     }
    //     if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    //       generator()->log() << "repairing remnant " << r1->id()
    // 			 << " and " << r2->id() << endl;
    //     }
    //     r1->setMomentum(r1->parents()[0]->momentum() - q1->momentum());
    //     r2->setMomentum(r2->parents()[0]->momentum() - q2->momentum());
    //     if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    //       generator()->log() << "r1->momentum() = " << r1->momentum() << endl
    // 			 << "r2->momentum() = " << r2->momentum() << endl;
    //     }
    // //     r1->boost(boostRest); 
    // //     r1->boost(boostNewF);
    // //     r2->boost(boostRest);
    // //     r2->boost(boostNewF);

    // //    cout << "1st jet" << endl;
    //     while(abs(q1->children()[0]->id()) < 22) {
    //       if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) {
    // 	generator()->log() << (q1->momentum() - q1->children()[0]->momentum() 
    // 			       - q1->children()[1]->momentum())/GeV << endl;
    //       }
    //       q1 = q1->children()[0];
    //     }
    //     //    cout << "2nd jet" << endl;
    //     while(abs(q2->children()[0]->id()) < 22) {
    //       if ( HERWIG_DEBUG_LEVEL >= HwDebug::full_Shower ) { 
    // 	generator()->log() << (q2->momentum() - q2->children()[0]->momentum() 
    // 			       - q2->children()[1]->momentum())/GeV << endl;
    //       }
    //       q2 = q2->children()[0];
    //     }
    //     



// comment block 2

    
    //      * to book some histograms...

    //     static SampleHistogram Qhist(-20.0, 1000.0, 20.0);
    //     static SampleHistogram pLong(0.0, 5000.0, 100.0);
    //     static SampleHistogram pt(0.0, 1000.0, 20.0);
    //     static SampleHistogram Energy(-500.0, 5000.0, 100.0);
    // //     static SampleHistogram Qhist(-2.0, 100.0, 2.0);
    // //     static SampleHistogram pLong(0.0, 500.0, 10.0);
    // //     static SampleHistogram pt(0.0, 100.0, 2.0);
    // //     static SampleHistogram Energy(-100.0, 500.0, 10.0);
    //     Qhist += -toBoost[0]->momentum().m()/GeV;
    //     Qhist += -toBoost[1]->momentum().m()/GeV;
    //     pLong += toBoost[0]->momentum().z()/GeV;
    //     pLong += toBoost[1]->momentum().z()/GeV;
    //     pt += toBoost[0]->momentum().perp()/GeV;
    //     pt += toBoost[1]->momentum().perp()/GeV;
    //     Energy += toBoost[0]->momentum().e()/GeV;
    //     Energy += toBoost[1]->momentum().e()/GeV;
    //     Qhist.printGnuplot("Qhist.dat");
    //     pLong.printGnuplot("pLong.dat");
    //     pt.printGnuplot("pperp.dat");
    //     Energy.printGnuplot("energy.dat");

    //     static SampleHistogram Q_u_p(-20.0, 1000.0, 20.0);  
    //     static SampleHistogram Q_d_p(-20.0, 1000.0, 20.0);  
    //     static SampleHistogram Q_ub_p(-20.0, 1000.0, 20.0); 
    //     static SampleHistogram Q_db_p(-20.0, 1000.0, 20.0); 
    //     static SampleHistogram Q_u_pb(-20.0, 1000.0, 20.0); 
    //     static SampleHistogram Q_d_pb(-20.0, 1000.0, 20.0); 
    //     static SampleHistogram Q_ub_pb(-20.0, 1000.0, 20.0);
    //     static SampleHistogram Q_db_pb(-20.0, 1000.0, 20.0);
			   						     
    //     static SampleHistogram pt_u_p(0.0, 1000.0, 20.0); 
    //     static SampleHistogram pt_d_p(0.0, 1000.0, 20.0); 
    //     static SampleHistogram pt_ub_p(0.0, 1000.0, 20.0);
    //     static SampleHistogram pt_db_p(0.0, 1000.0, 20.0);
    //     static SampleHistogram pt_u_pb(0.0, 1000.0, 20.0);
    //     static SampleHistogram pt_d_pb(0.0, 1000.0, 20.0);
    //     static SampleHistogram pt_ub_pb(0.0, 1000.0, 20.0);
    //     static SampleHistogram pt_db_pb(0.0, 1000.0, 20.0);
			   						     
    //     static SampleHistogram pl_u_p(0.0, 5000.0, 100.0); 
    //     static SampleHistogram pl_d_p(0.0, 5000.0, 100.0); 
    //     static SampleHistogram pl_ub_p(0.0, 5000.0, 100.0);
    //     static SampleHistogram pl_db_p(0.0, 5000.0, 100.0);
    //     static SampleHistogram pl_u_pb(0.0, 5000.0, 100.0);
    //     static SampleHistogram pl_d_pb(0.0, 5000.0, 100.0);
    //     static SampleHistogram pl_ub_pb(0.0, 5000.0, 100.0);
    //     static SampleHistogram pl_db_pb(0.0, 5000.0, 100.0);
			   						     
    //     static SampleHistogram e_u_p(-500.0, 5000.0, 100.0);  
    //     static SampleHistogram e_d_p(-500.0, 5000.0, 100.0);  
    //     static SampleHistogram e_ub_p(-500.0, 5000.0, 100.0); 
    //     static SampleHistogram e_db_p(-500.0, 5000.0, 100.0); 
    //     static SampleHistogram e_u_pb(-500.0, 5000.0, 100.0); 
    //     static SampleHistogram e_d_pb(-500.0, 5000.0, 100.0); 
    //     static SampleHistogram e_ub_pb(-500.0, 5000.0, 100.0);
    //     static SampleHistogram e_db_pb(-500.0, 5000.0, 100.0);

    // //     static SampleHistogram Q_u_p(-2.0, 100.0, 2.0);  
    // //     static SampleHistogram Q_d_p(-2.0, 100.0, 2.0);  
    // //     static SampleHistogram Q_ub_p(-2.0, 100.0, 2.0); 
    // //     static SampleHistogram Q_db_p(-2.0, 100.0, 2.0); 
    // //     static SampleHistogram Q_u_pb(-2.0, 100.0, 2.0); 
    // //     static SampleHistogram Q_d_pb(-2.0, 100.0, 2.0); 
    // //     static SampleHistogram Q_ub_pb(-2.0, 100.0, 2.0);
    // //     static SampleHistogram Q_db_pb(-2.0, 100.0, 2.0);
			   						     
    // //     static SampleHistogram pt_u_p(0.0, 100.0, 2.0); 
    // //     static SampleHistogram pt_d_p(0.0, 100.0, 2.0); 
    // //     static SampleHistogram pt_ub_p(0.0, 100.0, 2.0);
    // //     static SampleHistogram pt_db_p(0.0, 100.0, 2.0);
    // //     static SampleHistogram pt_u_pb(0.0, 100.0, 2.0);
    // //     static SampleHistogram pt_d_pb(0.0, 100.0, 2.0);
    // //     static SampleHistogram pt_ub_pb(0.0, 100.0, 2.0);
    // //     static SampleHistogram pt_db_pb(0.0, 100.0, 2.0);
			   						     
    // //     static SampleHistogram pl_u_p(0.0, 500.0, 10.0); 
    // //     static SampleHistogram pl_d_p(0.0, 500.0, 10.0); 
    // //     static SampleHistogram pl_ub_p(0.0, 500.0, 10.0);
    // //     static SampleHistogram pl_db_p(0.0, 500.0, 10.0);
    // //     static SampleHistogram pl_u_pb(0.0, 500.0, 10.0);
    // //     static SampleHistogram pl_d_pb(0.0, 500.0, 10.0);
    // //     static SampleHistogram pl_ub_pb(0.0, 500.0, 10.0);
    // //     static SampleHistogram pl_db_pb(0.0, 500.0, 10.0);
			   						     
    // //     static SampleHistogram e_u_p(-100.0, 500.0, 10.0);  
    // //     static SampleHistogram e_d_p(-100.0, 500.0, 10.0);  
    // //     static SampleHistogram e_ub_p(-100.0, 500.0, 10.0); 
    // //     static SampleHistogram e_db_p(-100.0, 500.0, 10.0); 
    // //     static SampleHistogram e_u_pb(-100.0, 500.0, 10.0); 
    // //     static SampleHistogram e_d_pb(-100.0, 500.0, 10.0); 
    // //     static SampleHistogram e_ub_pb(-100.0, 500.0, 10.0);
    // //     static SampleHistogram e_db_pb(-100.0, 500.0, 10.0);

    //     if (toBoost[ia]->id() == 2) {
    //       Q_u_p += -toBoost[ia]->momentum().m()/GeV;
    //       pl_u_p += toBoost[ia]->momentum().z()/GeV;
    //       pt_u_p += toBoost[ia]->momentum().perp()/GeV;
    //       e_u_p += toBoost[ia]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ia]->id() == 1) {
    //       Q_d_p += -toBoost[ia]->momentum().m()/GeV;
    //       pl_d_p += toBoost[ia]->momentum().z()/GeV;
    //       pt_d_p += toBoost[ia]->momentum().perp()/GeV;
    //       e_d_p += toBoost[ia]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ia]->id() == -2) {
    //       Q_ub_p += -toBoost[ia]->momentum().m()/GeV;
    //       pl_ub_p += toBoost[ia]->momentum().z()/GeV;
    //       pt_ub_p += toBoost[ia]->momentum().perp()/GeV;
    //       e_ub_p += toBoost[ia]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ia]->id() == -1) {
    //       Q_db_p += -toBoost[ia]->momentum().m()/GeV;
    //       pl_db_p += toBoost[ia]->momentum().z()/GeV;
    //       pt_db_p += toBoost[ia]->momentum().perp()/GeV;
    //       e_db_p += toBoost[ia]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ib]->id() == 2) {
    //       Q_u_pb += -toBoost[ib]->momentum().m()/GeV;
    //       pl_u_pb += toBoost[ib]->momentum().z()/GeV;
    //       pt_u_pb += toBoost[ib]->momentum().perp()/GeV;
    //       e_u_pb += toBoost[ib]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ib]->id() == 1) {
    //       Q_d_pb += -toBoost[ib]->momentum().m()/GeV;
    //       pl_d_pb += toBoost[ib]->momentum().z()/GeV;
    //       pt_d_pb += toBoost[ib]->momentum().perp()/GeV;
    //       e_d_pb += toBoost[ib]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ib]->id() == -2) {
    //       Q_ub_pb += -toBoost[ib]->momentum().m()/GeV;
    //       pl_ub_pb += toBoost[ib]->momentum().z()/GeV;
    //       pt_ub_pb += toBoost[ib]->momentum().perp()/GeV;
    //       e_ub_pb += toBoost[ib]->momentum().e()/GeV;
    //     }
    //     if (toBoost[ib]->id() == -1) {
    //       Q_db_pb += -toBoost[ib]->momentum().m()/GeV;
    //       pl_db_pb += toBoost[ib]->momentum().z()/GeV;
    //       pt_db_pb += toBoost[ib]->momentum().perp()/GeV;
    //       e_db_pb += toBoost[ib]->momentum().e()/GeV;
    //     }

    //     Q_u_p.printGnuplot("fl-0.dat");  
    //     Q_d_p.printGnuplot("fl-1.dat");  
    //     Q_ub_p.printGnuplot("fl-2.dat"); 
    //     Q_db_p.printGnuplot("fl-3.dat"); 
    //     Q_u_pb.printGnuplot("fl-4.dat"); 
    //     Q_d_pb.printGnuplot("fl-5.dat"); 
    //     Q_ub_pb.printGnuplot("fl-6.dat");
    //     Q_db_pb.printGnuplot("fl-7.dat");
    //     pt_u_p.printGnuplot("fl-8.dat"); 
    //     pt_d_p.printGnuplot("fl-9.dat"); 
    //     pt_ub_p.printGnuplot("fl-10.dat");
    //     pt_db_p.printGnuplot("fl-11.dat");
    //     pt_u_pb.printGnuplot("fl-12.dat");
    //     pt_d_pb.printGnuplot("fl-13.dat");
    //     pt_ub_pb.printGnuplot("fl-14.dat");
    //     pt_db_pb.printGnuplot("fl-15.dat");				     
    //     pl_u_p.printGnuplot("fl-16.dat"); 
    //     pl_d_p.printGnuplot("fl-17.dat"); 
    //     pl_ub_p.printGnuplot("fl-18.dat");
    //     pl_db_p.printGnuplot("fl-19.dat");
    //     pl_u_pb.printGnuplot("fl-20.dat");
    //     pl_d_pb.printGnuplot("fl-21.dat");
    //     pl_ub_pb.printGnuplot("fl-22.dat");
    //     pl_db_pb.printGnuplot("fl-23.dat");
    //     e_u_p.printGnuplot("fl-24.dat");  
    //     e_d_p.printGnuplot("fl-25.dat");  
    //     e_ub_p.printGnuplot("fl-26.dat"); 
    //     e_db_p.printGnuplot("fl-27.dat"); 
    //     e_u_pb.printGnuplot("fl-28.dat"); 
    //     e_d_pb.printGnuplot("fl-29.dat"); 
    //     e_ub_pb.printGnuplot("fl-30.dat");
    //     e_db_pb.printGnuplot("fl-31.dat");
    //    


// comment block 3

    ///////////////////////////////////////////////////////////////////////////////

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
