// -*- C++ -*-
//
#ifndef HERWIG_MyKinematics_H
#define HERWIG_MyKinematics_H
//
// This is the declaration of the <!id>Kinematics<!!id> class.
//
// CLASSDOC SUBSECTION Description:
// 
// This is a pure static class which provides some useful methods
// for kinematics computation, as the two body decays.
//
// NB) Other useful kinematical methods (and probably even those
//     implemented in Kinematics class!) are in
//         Pythia7/Utilities/UtilityBase.h
//

#include "Herwig++/Config/Herwig.h"
#include "Pythia7/CLHEPWrap/ThreeVector.h"
#include "Pythia7/CLHEPWrap/LorentzRotation.h"
#include "Pythia7/Repository/EventGenerator.h"


namespace Herwig {


  using namespace Pythia7;


  class MyKinematics {

  public:

    static void twoBodyDecay(const Lorentz5Momentum & p, 
			     const double m1, const double m2,
			     const Vector3 & unitDir1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);
    static void twoBodyDecay(const Lorentz5Momentum & p, 
			     const double m1, const double m2,
			     const double cosThetaStar1, const double phiStar1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);
    // Given in input: the 5-momentum (p) of the decay particle in the LAB frame;
    // the masses of the two decay products (m1,m2); either the direction of the
    // first of the two decay products (unitDir1), or its cosine of the polar angle
    // (cosThetaStar1) and the azimuthal angle (phiStar1), in the CM system of the
    // decaying particle; the methods gives the 5-momenta of the two decay products 
    // in the LAB frame (p1,p2).

    static inline Energy CMMomentum(const Energy M, 
				    const Energy m1, 
				    const Energy m2);
    // For the two body decay  M -> m1 + m2  it gives the module of the 3-momentum
    // of the decay product in the rest frame of M.

    static inline Vector3 unitDirection(const double cosTheta, const double phi);
    // It returns the unit 3-vector with the given  cosTheta  and  phi. 
    
    static LorentzRotation rotation(const Lorentz5Momentum, const double, const double);
    // Create matrix that rotates 5 mom to z axis, then rotates by phi 
    // (where phi is given as two double args, cos-phi and sin-phi

    static void generateAngles(tEGPtr, double &, double &);
    // A routine to randomly generate two angles. The first double argument is
    // a cos of an angle (between -1 and 1) and the second double is an
    // angle between 0 and 2 pi
  private:

    MyKinematics();
    MyKinematics(const MyKinematics & x);
    MyKinematics & operator=(const MyKinematics & x);
    // pure static class

  };


}

#include "Kinematics.icc"

#endif /* HERWIG_Kinematics_H */

