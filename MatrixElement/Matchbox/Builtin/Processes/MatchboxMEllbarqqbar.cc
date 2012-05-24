// -*- C++ -*-
//
// MatchboxMEllbarqqbar.cc is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxMEllbarqqbar.h"

using namespace Herwig;
using namespace SpinorHelicity;

void MatchboxMEllbarqqbar::doinit(const HandlerBase& h) {

  // make sure this is initialized, as it may not have happened
  // before we get here due to the preinitRegister fiddling
  // for the dipoles
  h.generator()->standardModel()->init();

  // Vff = ie \gamma^\mu ( v_V - a_V\gamma^5 ) with e > 0
  // v_Z = (I_3 - 2 s^2 Q)/(2 s c) , a_Z = I_3/(2 s c)
  // v_\gamma = -Q , a_\gamma = 0

  BMass = h.getParticleData(ParticleID::Z0)->mass();
  BWidth = h.getParticleData(ParticleID::Z0)->width();

  Nc = h.generator()->standardModel()->Nc();
  alphaS = h.generator()->standardModel()->alphaS();

  double s2tW = h.generator()->standardModel()->sin2ThetaW();
  double thetaWFactor = 4.*sqrt(s2tW*(1.-s2tW)); 
  double alphaEM = h.generator()->standardModel()->alphaEM();
  double alphaEMFactor = sqrt(4.*Constants::pi*alphaEM);

  el = - h.generator()->standardModel()->ee() * alphaEMFactor;
  vl = h.generator()->standardModel()->ve()  * alphaEMFactor/thetaWFactor;
  al = h.generator()->standardModel()->ae() * alphaEMFactor/thetaWFactor;

  quarkQEDCouplings = make_pair(- h.generator()->standardModel()->eu() * alphaEMFactor,
				- h.generator()->standardModel()->ed() * alphaEMFactor);
  quarkVectorCouplings = make_pair(h.generator()->standardModel()->vu() * alphaEMFactor/thetaWFactor,
				   h.generator()->standardModel()->vd() * alphaEMFactor/thetaWFactor);
  quarkAxialCouplings = make_pair(h.generator()->standardModel()->au() * alphaEMFactor/thetaWFactor,
				  h.generator()->standardModel()->ad() * alphaEMFactor/thetaWFactor);

  nPoints(4);

}

void MatchboxMEllbarqqbar::persistentOutput(PersistentOStream & os) const {

  os << ounit(BMass,GeV) << ounit(BWidth,GeV) 
     << Nc << alphaS
     << el << vl << al 
     << quarkQEDCouplings << quarkVectorCouplings << quarkAxialCouplings;

}
						      
void MatchboxMEllbarqqbar::persistentInput(PersistentIStream & is) {

  is >> iunit(BMass,GeV) >> iunit(BWidth,GeV) 
     >> Nc >> alphaS
     >> el >> vl >> al 
     >> quarkQEDCouplings >> quarkVectorCouplings >> quarkAxialCouplings;

  nPoints(4);

}

void MatchboxMEllbarqqbar::prepare(const Lorentz5Momentum& pl, const Lorentz5Momentum& plbar,
				   const Lorentz5Momentum& pq, const Lorentz5Momentum& pqbar,
				   Energy2 sHat, 
				   cPDPtr lData, cPDPtr lbarData, 
				   cPDPtr qData, cPDPtr qbarData) const {

  amplitudeScale(sqrt(sHat));

  momentum(0,pl,false,lData->mass());
  momentum(1,plbar,false,lbarData->mass());

  momentum(2,pq,false,qData->mass()); 
  momentum(3,pqbar,false,qbarData->mass());

  eq =
    abs(qData->id()) % 2 == 0 ?
    quarkQEDCouplings.first :
    quarkQEDCouplings.second;

  vq =
    abs(qData->id()) % 2 == 0 ?
    quarkVectorCouplings.first :
    quarkVectorCouplings.second;

  aq =
    abs(qData->id()) % 2 == 0 ?
    quarkAxialCouplings.first :
    quarkAxialCouplings.second;

}

double MatchboxMEllbarqqbar::evaluateME2(bool photon) const {

  double xQ2 = (momentum(0)+momentum(1)).m2();
  double M2 = sqr(BMass/amplitudeScale());
  double G2 = sqr(BWidth/amplitudeScale());

  double BB =
    8*(4*al*aq*vl*vq*(-(invariant(0,3)*invariant(1,2)) + invariant(0,2)*invariant(1,3)) + 
     (invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) + 2*invariant(2,3)*mass(0)*mass(1) - 
        2*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3))*sqr(aq)*sqr(vl) + 
     (invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) + 2*invariant(2,3)*mass(0)*mass(1) + 
        2*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3))*sqr(vl)*sqr(vq) + 
     sqr(al)*((invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) - 2*invariant(2,3)*mass(0)*mass(1) - 
           2*(invariant(0,1) - 4*mass(0)*mass(1))*mass(2)*mass(3))*sqr(aq) + 
        (invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) + 2*invariant(0,1)*mass(2)*mass(3) - 
	 2*mass(0)*mass(1)*(invariant(2,3) + 4*mass(2)*mass(3)))*sqr(vq)));

  double res =
    BB / ( sqr(xQ2-M2) + M2*G2 );

  if ( photon ) {

    double PP =
      8*(invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) + 2*invariant(2,3)*mass(0)*mass(1) + 
	 2*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3))*sqr(el)*sqr(eq);

    double PB = 
      8*el*eq*(al*aq*(-(invariant(0,3)*invariant(1,2)) + invariant(0,2)*invariant(1,3)) + 
	       vl*vq*(invariant(0,3)*invariant(1,2) + invariant(0,2)*invariant(1,3) + 2*invariant(2,3)*mass(0)*mass(1) + 
		      2*(invariant(0,1) + 4*mass(0)*mass(1))*mass(2)*mass(3)));

    res += PP / sqr(xQ2) + 2.*(1.-M2/xQ2)*PB / ( sqr(xQ2-M2) + M2*G2 );

  }

  return Nc*res;

}
