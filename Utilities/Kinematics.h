// -*- C++ -*-
//
#ifndef HERWIG_Kinematics_H
#define HERWIG_Kinematics_H
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
//         ThePEG/Utilities/UtilityBase.h
//

#include "Herwig++/Config/Herwig.h"
#include "ThePEG/CLHEPWrap/ThreeVector.h"
#include "ThePEG/CLHEPWrap/LorentzRotation.h"

namespace Herwig {


  using namespace ThePEG;


  class Kinematics {

  public:

    static void twoBodyDecay(const Lorentz5Momentum & p, 
			     const double m1, const double m2,
			     const Vector3 & unitDir1,
			     Lorentz5Momentum & p1, Lorentz5Momentum & p2);
    static void twoBodyDecay(const Lorentz5Momentum &, 
			     const double, const double,
			     const double, const double,
			     Lorentz5Momentum &, Lorentz5Momentum &);
    /***
     * Given in input: the 5-momentum (p) of the decay particle in the LAB 
     * frame; the masses of the two decay products (m1,m2); either the 
     * direction of the first of the two decay products (unitDir1), or its 
     * cosine of the polar angle (cosThetaStar1) and the azimuthal angle 
     * (phiStar1), in the CM system of the decaying particle; the methods 
     * gives the 5-momenta of the two decay products in the LAB frame (p1,p2).
     ***/
    static void threeBodyDecay(Lorentz5Momentum p0, Lorentz5Momentum &p1, 
			       Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			       double (*fcn)(double*) = NULL);
    /***
     * As the name implies, this takes the momentum p0 and does a flat three
     * body decay into p1..p3. The argument fcn is used to add additional
     * weights. If it is not used, the default is just flat in phasespace.
     ***/
    static void fourBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			      Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			      Lorentz5Momentum &p4);
    /***
     * Again, as the name implies, this is an isotropic four-body decay
     ***/
    static void fiveBodyDecay(Lorentz5Momentum  p0, Lorentz5Momentum &p1,
			      Lorentz5Momentum &p2, Lorentz5Momentum &p3,
			      Lorentz5Momentum &p4, Lorentz5Momentum &p5);
    /***
     * Again as the name implies, this is an isotropic five-body decay 
     ***/

    static inline Energy pstarTwoBodyDecay(const Energy M, 
					   const Energy m1, const Energy m2);
    /***
     * For the two body decay  M -> m1 + m2  it gives the module of the 
     * 3-momentum of the decay product in the rest frame of M.
     ***/

    static inline Vector3 unitDirection(const double, const double);
    // It returns the unit 3-vector with the given  cosTheta  and  phi. 

    static Energy CMMomentum(const Energy M, 
				  const Energy m1, 
				  const Energy m2);
    // This returns the CMMomentum of a two body decay, given M, m1, m2

    static LorentzRotation rotation(const Lorentz5Momentum p,
				    const double cp, 
				    const double sp);
    // This returns the Matrix that rotates along the direction p by
    // phi (cp = cos(phi), sp = sin(phi))

    static void generateAngles(double &, double &);
    // This just generates angles. First flat -1..1, second flat 0..2Pi

  private:

    Kinematics();
    Kinematics(const Kinematics & x);
    Kinematics & operator=(const Kinematics & x);
    // pure static class

  };


}

#include "Kinematics.icc"

#endif /* HERWIG_Kinematics_H */

