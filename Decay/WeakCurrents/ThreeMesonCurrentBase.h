// -*- C++ -*-
//
// ThreeMesonCurrentBase.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ThreeMesonCurrentBase_H
#define HERWIG_ThreeMesonCurrentBase_H
// This is the declaration of the ThreeMesonCurrentBase class.

#include "WeakDecayCurrent.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  This is the base class for the three meson decays of the weak current.
 *  It is designed so that the currents for the following modes can be implemented
 *  in classes inheriting from this
 * - \f$    \pi^-  \pi^-    \pi^+ \f$, (imode=0)
 * - \f$    \pi^0  \pi^0    \pi^- \f$, (imode=1)
 * - \f$    K^-   \pi^-    K^+ \f$, (imode=2)
 * - \f$    K^0   \pi^-    \bar{K}^0\f$, (imode=3)
 * - \f$    K^-   \pi^0    K^0 \f$, (imode=4)
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=5)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=6)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=7)
 * - \f$    \pi^-  \pi^0    \eta \f$, (imode=8)
 *
 * obviously there are other modes with three pseudoscalar mesons for the decay
 * of the weak current but this model original came from \f$\tau\f$ decay where
 * these are the only modes. However one case which is important is the inclusion
 * of the mixing in the neutral kaon sector for which we include the additional
 * currents
 * - \f$    K^0_S \pi^- K^0_S\f$, (imode=9)
 * - \f$    K^0_L \pi^- K^0_L\f$, (imode=10)
 * - \f$    K^0_S \pi^- K^0_L\f$, (imode=11)
 *
 *  In this case the current is given by
 *  \f[ J^\mu = \left(g^{\mu\nu}-\frac{q^\mu q^\nu}{q^2}\right)
 *   \left[F_1(p_2-p_3)^\mu +F_2(p_3-p_1)^\mu+F_3(p_1-p_2)^\mu\right]
 *  +q^\mu F_4
 *  +F_5\epsilon^{\mu\alpha\beta\gamma}p_1^\alpha p_2^\beta p_3^\gamma
 *  \f]
 * where
 * - \f$p_{1,2,3}\f$ are the momenta of the mesons in the order given above.
 * - \f$F_1,F_2,F_3,F_4,F_5\f$ are the form factors which must be 
 *  calculated in the calculateFormFactors member which should be implemented
 * in classes inheriting from this.
 *
 * @see WeakDecayCurrent.
 *  
 * \author Peter Richardson
 *
 */
class ThreeMesonCurrentBase: public WeakDecayCurrent {

public:

  /**
   * Default constructor
   */
  ThreeMesonCurrentBase();

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

public:


  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(const int imode,const int ichan,Energy & scale,
	  const ParticleVector & decay,DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Checks the mesons against the list.
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * The particles produced by the current. This returns the mesons for the mode.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

protected:

  /**
   * can a particular decayer handle this type of mode
   * @param imode The mode number as given above
   * @return Whether this mode can be handled.
   */
  virtual bool acceptMode(int imode) const=0;

  /**
   * Helper class for form factors
   */
  struct FormFactors {

    /**
     * @param F1 The \f$F_1\f$ form factor
     */
    complex<InvEnergy>  F1;
    
    /**
     * @param F2 The \f$F_2\f$ form factor
     */
    complex<InvEnergy>  F2;
    
    /**
     * @param F3 The \f$F_3\f$ form factor
     */
    complex<InvEnergy>  F3; 
    
    /**
     * @param F4 The \f$F_4\f$ form factor
     */
    complex<InvEnergy>  F4;
    
    /**
     * @param F5 The \f$F_5\f$ form factor
     */
    complex<InvEnergy3> F5;

    /**
     *  Constructor
     * @param f1 The \f$F_1\f$ form factor
     * @param f2 The \f$F_2\f$ form factor
     * @param f3 The \f$F_3\f$ form factor
     * @param f4 The \f$F_4\f$ form factor
     * @param f5 The \f$F_5\f$ form factor
     */    
    FormFactors(complex<InvEnergy>  f1 = InvEnergy(), 
		complex<InvEnergy>  f2 = InvEnergy(),
		complex<InvEnergy>  f3 = InvEnergy(),
		complex<InvEnergy>  f4 = InvEnergy(),
		complex<InvEnergy3> f5 = InvEnergy3())
      : F1(f1), F2(f2), F3(f3), F4(f4), F5(f5) {}
  };

  /**
   * Calculate the form factor for the current.
   * @param ichan The phase space channel
   * @param imode The mode
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   */
  virtual FormFactors calculateFormFactors(const int ichan, const int imode,
					   Energy2 q2,
					   Energy2 s1, Energy2 s2, Energy2 s3) const = 0;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ThreeMesonCurrentBase & operator=(const ThreeMesonCurrentBase &);

};

}

#endif /* HERWIG_ThreeMesonCurrentBase_H */
