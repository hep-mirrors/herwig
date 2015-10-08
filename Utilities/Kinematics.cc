// -*- C++ -*-
//
// Kinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Kinematics class.
//

#include "Kinematics.h"
#include <ThePEG/Vectors/Lorentz5Vector.h>
#include <ThePEG/Vectors/LorentzVector.h>
#include <ThePEG/Vectors/LorentzRotation.h>
#include <ThePEG/Repository/EventGenerator.h>
#include <ThePEG/Repository/CurrentGenerator.h>
#include <ThePEG/EventRecord/Event.h>

using namespace Herwig;
using namespace ThePEG;

bool Kinematics::twoBodyDecay(const Lorentz5Momentum & p, 
           const Energy m1, const Energy m2,
           const Axis & unitDir1,
           Lorentz5Momentum & p1, Lorentz5Momentum & p2) {
      Energy min=p.mass();
      if ( min >= m1 + m2  &&  m1 >= ZERO  &&  m2 >= ZERO  ) {
  Momentum3 pstarVector = unitDir1 * pstarTwoBodyDecay(min,m1,m2);
  p1 = Lorentz5Momentum(m1, pstarVector);
  p2 = Lorentz5Momentum(m2,-pstarVector);
  // boost from CM to LAB
  Boost bv = p.boostVector();
  double gammarest = p.e()/p.mass();
  p1.boost( bv, gammarest );
  p2.boost( bv, gammarest );
  return true;
      }
      return false;
    }

/*****
 * This function, as the name implies, performs a three body decay. The decay
 * products are distributed uniformly in all three directions.
 ****/
bool Kinematics::threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
				Lorentz5Momentum &p2, Lorentz5Momentum &p3,
				double (*fcn)(Energy2,Energy2,Energy2,InvEnergy4)) {
  // Variables needed in calculation...named same as fortran version
  Energy a = p0.mass() + p1.mass();
  Energy b = p0.mass() - p1.mass();
  Energy c = p2.mass() + p3.mass();
  
  if(b < c) {
     CurrentGenerator::log() 
       << "Kinematics::threeBodyDecay() phase space problem\n"
       << p0.mass()/GeV << " -> "
       << p1.mass()/GeV << ' '
       << p2.mass()/GeV << ' '
       << p3.mass()/GeV << '\n';
     return false;
  }
  
  Energy d = abs(p2.mass()-p3.mass());
  Energy2 aa = sqr(a); 
  Energy2 bb = sqr(b); 
  Energy2 cc = sqr(c); 
  Energy2 dd = sqr(d); 
  Energy2 ee = (b-c)*(a-d);
  
  Energy2 a1 = 0.5 * (aa+bb);
  Energy2 b1 = 0.5 * (cc+dd);
  InvEnergy4 c1 = 4./(sqr(a1-b1));
  
  Energy2 ff; 
  double ww; 
  Energy4 pp,qq,rr;
  // Choose mass of subsystem 23 with prescribed distribution
  const unsigned int MAXTRY = 100;
  unsigned int ntry=0;
  do {
    // ff is the mass squared of the 23 subsystem
    ff = UseRandom::rnd()*(cc-bb)+bb;
    
    // pp is ((m0+m1)^2 - m23^2)((m0-m1)^2-m23)
    pp = (aa-ff)*(bb-ff);
    
    // qq is ((m2+m3)^2 - m23^2)(|m2-m3|^2-m23^2)
    qq = (cc-ff)*(dd-ff);
    
    // weight
    ww = (fcn != NULL) ? (*fcn)(ff,a1,b1,c1) : 1.0;
    ww = sqr(ww);
    rr = ee*ff*UseRandom::rnd();
    ++ntry;
  } 
  while(pp*qq*ww < rr*rr && ntry < MAXTRY );
  if(ntry >= MAXTRY) {
    CurrentGenerator::log() << "Kinematics::threeBodyDecay can't generate momenta" 
			    << " after " << MAXTRY << " attempts\n";
    return false;
  }

  // ff is the mass squared of subsystem 23
  // do 2 body decays 0->1+23, 23->2+3
  double CosAngle, AzmAngle;
  Lorentz5Momentum p23;
  
  p23.setMass(sqrt(ff));
  
  generateAngles(CosAngle,AzmAngle);
  bool status = twoBodyDecay(p0,p1.mass(),p23.mass(),CosAngle,AzmAngle,p1,p23);
  
  generateAngles(CosAngle,AzmAngle);
  status &= twoBodyDecay(p23,p2.mass(),p3.mass(),CosAngle,AzmAngle,p2,p3);
  return status;
}
