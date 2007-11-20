// -*- C++ -*-
//
// Kinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
#ifndef HERWIG_Kinematics_H
#define HERWIG_Kinematics_H

// This is the declaration of the Kinematics class.

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/Vectors/LorentzRotation.h"

namespace Herwig {

  using namespace ThePEG;

  /** \ingroup Utilities
   *  This is a pure static class which provides some useful methods
   *  for kinematics computation, as the two body decays.
   * 
   *  NB) For other useful kinematical methods (and probably even those
   *      implemented in Kinematics class!):
   *          @see UtilityBase
   */
  class Kinematics {

  public:

    /**
     *  Calculate the momenta for a two body decay
     * The return value indicates success or failure.
     * @param p The momentum of the decaying particle
     * @param m1 The mass of the first decay product
     * @param m2 The mass of the seocnd decay product
     * @param unitDir1 Direction for the products in the rest frame of
     * the decaying particle
     * @param p1 The momentum of the first decay product
     * @param p2 The momentum of the second decay product
     */
    static bool twoBodyDecay(const Lorentz5Momentum & p, 
			     const Energy m1, const Energy m2,
			     const Axis & unitDir1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);

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
    static bool twoBodyDecay(const Lorentz5Momentum & p, 
			     const Energy m1, const Energy m2,
			     const double cosThetaStar1, 
			     const double phiStar1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);

    /**
     * As the name implies, this takes the momentum p0 and does a flat three
     * body decay into p1..p3. The argument fcn is used to add additional
     * weights. If it is not used, the default is just flat in phasespace.
     * The return value indicates success or failure.
     */
    static bool threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       double (*fcn)(Energy2,Energy2,Energy2,InvEnergy4) = NULL);


    /**
     * Again as the name implies, this is an isotropic four-body decay.
     * The return value indicates success or failure.
     */
    static bool fourBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			      Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			      Lorentz5Momentum &p4);

    /**
     * Again as the name implies, this is an isotropic five-body decay.
     * The return value indicates success or failure.
     */
    static bool fiveBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			      Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			      Lorentz5Momentum &p4, Lorentz5Momentum &p5);

    /**
     * For the two body decay  M -> m1 + m2  it gives the module of the 
     * 3-momentum of the decay product in the rest frame of M.
     */
    static inline Energy pstarTwoBodyDecay(const Energy M, 
					   const Energy m1, const Energy m2);

    /**
     * It returns the unit 3-vector with the given  cosTheta  and  phi.
     */
    static inline Axis unitDirection(const double, const double);

    /**
     * This returns the CMMomentum of a two body decay, given M, m1, m2.
     */
    static Energy CMMomentum(const Energy M, 
			     const Energy m1, 
			     const Energy m2);

    /**
     * This just generates angles. First flat -1..1, second flat 0..2Pi
     */
    static void generateAngles(double &, double &);

  private:

    /**
     * Pure static class so default constructor private
     */
    Kinematics();

    /**
     * Pure static class so copy constructor private
     */
    Kinematics(const Kinematics & x);

    /**
     * The assignment operator is private and must never be called.
     * In fact, it should not even be implemented.
     */
    Kinematics & operator=(const Kinematics & x);

  };

}

#include "Kinematics.icc"

#endif /* HERWIG_Kinematics_H */

