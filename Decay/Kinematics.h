// -*- C++ -*-
#ifndef HERWIG_MyKinematics_H
#define HERWIG_MyKinematics_H
//
// This is the declaration of the Kinematics class.
//
#include "Herwig++/Config/Herwig.h"
#include <ThePEG/CLHEPWrap/ThreeVector.h>
#include <ThePEG/CLHEPWrap/LorentzRotation.h>
#include <ThePEG/Repository/EventGenerator.h>


namespace Herwig {
  using namespace ThePEG;

/** \ingroup Decay
 *
 * This is a pure static class which provides some useful methods
 * for kinematics computation, such as the two body decays.
 *
 * NB) Other useful kinematical methods (and probably even those
 *     implemented in Kinematics class!) are in
 *         ThePEG/Utilities/UtilityBase.h
 *
 * @see UtilityBase
 */

class MyKinematics {

  public:

  /**
   * Generate a two-body decay.
   * @param p The input momentum of the decaying particle
   * @param m1 The mass of the first decay product
   * @param m2 The mass of the second decay product
   * @param unitDir1 The direction of the first decay product.
   * @param p1 The output momentum of the first decay product.
   * @param p2 The output momentum of the second decay product.
   */
  static void twoBodyDecay(const Lorentz5Momentum & p, 
			   const double m1, const double m2,
			   const Vector3 & unitDir1,
			   Lorentz5Momentum & p1, Lorentz5Momentum & p2);
  
  /**
   * Generate a two-body decya
   * @param p The input momentum of the decaying particle
   * @param m1 The mass of the first decay product
   * @param m2 The mass of the second decay product
   * @param cosThetaStar1 the cosine of the polar angle for the first decay product.
   * @param phiStar1 The azimuthal angle of the first decay product.
   * @param p1 The output momentum of the first decay product.
   * @param p2 The output momentum of the second decay product.
   */
    static void twoBodyDecay(const Lorentz5Momentum & p, 
			     const double m1, const double m2,
			     const double cosThetaStar1, const double phiStar1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);

    /**
     * Momentum of the decay products in a two body decay
     * @param M mass of the decaying particle.
     * @param m1 The mass of the first  decay product.
     * @param m2 The mass of the second decay product.
     * @return The momentum of the decay products in the rest-frame of the decay
     *         particle.
     */
    static inline Energy CMMomentum(const Energy M, 
				    const Energy m1, 
				    const Energy m2);

    /**
     * It returns the unit 3-vector with the given  cosTheta  and  phi.
     * @param cosTheta The cosine of the polar angle.
     * @param phi The azimuthal angle
     * @return The unit vector. 
     */
    static inline Vector3 unitDirection(const double cosTheta, const double phi);
    
    /**
     * Create matrix that rotates a 5-momentum  to z-axis, then rotates by \f$\phi\f$ 
     * (where \f$\phi\f$ is given as two double args, \f$\cos\phi\f$ and \f$\sin\phi\f$.
     * @param mom The 5-momentum
     * @param cos Cosine of the angle  \f$\phi\f$, \f$\cos\phi\f$.
     * @param sin Sine of the angle  \f$\phi\f$, \f$\sin\phi\f$.
     */
    static LorentzRotation rotation(const Lorentz5Momentum mom, const double cos,
				    const double sin);

    /**
     * A routine to randomly generate two angles.
     * @param eg Pointer to an event generator to allow access to the random
     *                   number generator.
     * @param cosTheta The output value of the cosine of the polar angle between [-1,1]
     * @param phi The azimuthal angle between 0 and \f$2\pi\f$.
     */
    static void generateAngles(tEGPtr eg, double & cosTheta, double & phi);

  private:

    /**
     * pure static class, can't be created.
     */
    MyKinematics();

    /**
     * pure static class, can't be copied.
     */
    MyKinematics(const MyKinematics & x);

    /**
     * pure static class, can't be assigned.
     */
    MyKinematics & operator=(const MyKinematics & x);

  };


}

#include "Kinematics.icc"

#endif /* HERWIG_Kinematics_H */

