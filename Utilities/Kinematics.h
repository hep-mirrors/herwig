// -*- C++ -*-
//
// Kinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef HERWIG_Kinematics_H
#define HERWIG_Kinematics_H

// This is the declaration of the Kinematics class.

#include <ThePEG/Config/ThePEG.h>
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Repository/UseRandom.h"
#include <ThePEG/Vectors/Lorentz5Vector.h>

namespace Herwig {

  using namespace ThePEG;

  /** \ingroup Utilities
   *  This is a namespace which provides some useful methods
   *  for kinematics computation, as the two body decays.
   * 
   *  NB) For other useful kinematical methods (and probably even those
   *      implemented in Kinematics class!):
   *          @see UtilityBase
   */
  namespace Kinematics {

    /**
     *  Calculate the momenta for a two body decay
     * The return value indicates success or failure.
     * @param p The momentum of the decaying particle
     * @param m1 The mass of the first decay product
     * @param m2 The mass of the second decay product
     * @param unitDir1 Direction for the products in the rest frame of
     * the decaying particle
     * @param p1 The momentum of the first decay product
     * @param p2 The momentum of the second decay product
     */
    bool twoBodyDecay(const Lorentz5Momentum & p, 
			     const Energy m1, const Energy m2,
			     const Axis & unitDir1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);
    
    /**
     * It returns the unit 3-vector with the given  cosTheta  and  phi.
     */
    inline Axis unitDirection(const double cosTheta, const double phi) {
      return ( fabs( cosTheta ) <= 1.0  ? 
         Axis( cos(phi)*sqrt(1.0-cosTheta*cosTheta) , 
         sin(phi)*sqrt(1.0-cosTheta*cosTheta) , cosTheta) : Axis() );
    }

    /**
     *  Calculate the momenta for a two body decay
     * The return value indicates success or failure.
     * @param p The momentum of the decaying particle
     * @param m1 The mass of the first decay product
     * @param m2 The mass of the second decay product
     * @param cosThetaStar1 Polar angle in rest frame 
     * @param phiStar1 Azimuthal angle in rest frame
     * @param p1 The momentum of the first decay product
     * @param p2 The momentum of the second decay product
     */
    inline bool twoBodyDecay(const Lorentz5Momentum & p, 
			     const Energy m1, const Energy m2,
			     const double cosThetaStar1, 
			     const double phiStar1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2) {
      return twoBodyDecay(p,m1,m2,unitDirection(cosThetaStar1,phiStar1),p1,p2); 
    }

    /**
     * As the name implies, this takes the momentum p0 and does a flat three
     * body decay into p1..p3. The argument fcn is used to add additional
     * weights. If it is not used, the default is just flat in phasespace.
     * The return value indicates success or failure.
     */
    bool threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       double (*fcn)(Energy2,Energy2,Energy2,InvEnergy4) = NULL);

    /**
     * For the two body decay  M -> m1 + m2  it gives the module of the 
     * 3-momentum of the decay product in the rest frame of M.
     */
    inline Energy pstarTwoBodyDecay(const Energy M, 
				    const Energy m1, const Energy m2) {
      return ( M > ZERO &&  m1 >=ZERO && m2 >= ZERO  && M > m1+m2 ?  
	       Energy(sqrt(( sqr(M) - sqr(m1+m2) )*( sqr(M) - sqr(m1-m2) )) 
		      / (2.0*M) ) : ZERO); 
    }
    
    /**
     * This just generates angles. First flat -1..1, second flat 0..2Pi
     */
    inline void generateAngles(double & ct, double & az) {
      ct = UseRandom::rnd()*2.0 - 1.0;  // Flat from -1..1
      az = UseRandom::rnd()*Constants::twopi;   
    }
  }

}

#endif /* HERWIG_Kinematics_H */

