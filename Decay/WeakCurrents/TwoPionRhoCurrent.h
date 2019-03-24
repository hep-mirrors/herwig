// -*- C++ -*-
//
// TwoPionRhoCurrent.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoPionRhoCurrent_H
#define HERWIG_TwoPionRhoCurrent_H
// This is the declaration of the TwoPionRhoCurrent class.

#include "WeakCurrent.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/ResonanceHelpers.h"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 *  Weak current for the production of two mesons via the \f$\rho\f$ or \f$K^*\f$
 *  resonances.
 *  These currents are taken from tau decays.
 *
 *  The current takes the form
 *
 *  \f[J^\mu = \frac{\sqrt{2}}{\sum_k\alpha_k}\left((p_1-p_2)^\mu-\frac{(p_1-p_2)\cdot q}{q^2}q^\mu))\right)
 *   \sum_k \alpha_k B_{R_k}(q^2)
 *  \f]
 *  where
 *  - \f$p_{1,2}\f$ are the momenta of the outgoing mesons,
 *  - \f$q=p_1+p_2\f$,
 *  - \f$B_{R_k}(q^2)\f$ is the Breit-Wigner distribution for the intermediate vector
 *    meson \f$R_k\f$.
 *  - \f$\alpha_k\f$ is the weight for the resonance.
 *
 *   The Breit-Wigner term is summed over the \f$\rho\f$ or \f$K^*\f$ resonances that
 *   can contribute to a given decay.
 *
 *  The models of either Kuhn and Santamaria (Z. Phys. C48, 445 (1990))
 *  or Gounaris and Sakurai Phys. Rev. Lett. 21, 244 (1968) are supported for the
 *  shape of the Breit-Wigner distribution. The mixing parameters
 *  are taken from Phys.Rev.D61:112002,2000 (CLEO) for the decay \f$\pi^\pm\pi^0\f$ and
 *  the CLEO version of TAUOLA for the \f$K\pi\f$ decays.
 *
 * @see WeakCurrent.
 * 
 *  \author Peter Richardson
 *
 */
class TwoPionRhoCurrent: public WeakCurrent {

public:

  /**
   * Default constructor
   */
  TwoPionRhoCurrent();

  /** @name Methods for the construction of the phase space integrator. */
  //@{
  /**
   * Complete the construction of the decay mode for integration.classes inheriting
   * from this one.
   * This method is purely virtual and must be implemented in the classes inheriting
   * from WeakCurrent.
   * @param icharge   The total charge of the outgoing particles in the current.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
   * @param imode     The mode in the current being asked for.
   * @param mode      The phase space mode for the integration
   * @param iloc      The location of the of the first particle from the current in
   *                  the list of outgoing particles.
   * @param ires      The location of the first intermediate for the current.
   * @param phase     The prototype phase space channel for the integration.
   * @param upp       The maximum possible mass the particles in the current are
   *                  allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge, tcPDPtr resonance,
			  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			  unsigned int imode,PhaseSpaceModePtr mode,
			  unsigned int iloc,int ires,
			  PhaseSpaceChannel phase, Energy upp );

  /**
   * The particles produced by the current. This just returns the two pseudoscalar
   * mesons and the photon.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual tPDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This method is purely virtual and must be implemented in
   * all classes inheriting from this one.
   * @param resonance If specified only include terms with this particle
   * @param Itotal    If specified the total isospin of the current
   * @param I3        If specified the thrid component of isospin
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVectorE> 
  current(tcPDPtr resonance,
	  IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
	  const int imode, const int ichan,Energy & scale,
	  const tPDVector & outgoing,
	  const vector<Lorentz5Momentum> & momenta,
	  DecayIntegrator::MEOption meopt) const;

  /**
   * Accept the decay. Checks the particles are the allowed mode.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  TwoPionRhoCurrent & operator=(const TwoPionRhoCurrent &) = delete;

private:

  /**
   * \f$p\f$-wave breit wigner for form-factors
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param imodel Which of the two models for the Breit-Wigner shape to use.
   * @param ires   Which of the different multiplets to use.
   * @return The value of the Breit-Wigner distribution.
   */
  Complex BreitWigner(Energy2 q2, unsigned int imodel,
		      unsigned int ires) const {
    // calculate the BW
    if(imodel==0) {
      return Resonance::BreitWignerPWave(q2,_mass[ires],_width[ires],
					 _massa[ires],_massb[ires]);
    }
    else if(imodel==1) {
      return Resonance::BreitWignerGS(q2,_mass[ires],_width[ires],
				      _massa[ires],_massb[ires],
				      _h0[ires],_dh[ires],_hres[ires]);
    }
    else
      assert(false);
  }

private:

  /**
   * Weights for the different \f$\rho\f$ resonances in the current, \f$\alpha_k\f$.
   */
  //@{
  /**
   *  The Complex weight used in the calculation
   */
  vector<Complex> _piwgt;

  /**
   *  The magnitude for input
   */
  vector<double> _pimag;

  /**
   *  The phase for input
   */
  vector<double> _piphase;
  //@}

  /**
   * Model to use for the \f$\rho\f$ propagator.
   */
  int _pimodel;

  /**
   * Option not to use the physical masses and widths for the \f$\rho\f$.
   */
  bool _rhoparameters;

  /**
   * The masses of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhomasses;

  /**
   * The widths of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhowidths;

  /**
   * Parameters for the Breit-Wigners
   */
  //@{
  /**
   * The masses of the resonances
   */
  vector<Energy> _mass;

  /**
   * The widths of the resonances
   */
  vector<Energy> _width;

  /**
   * Masses of the decay products for the momentum calculation.
   */
  vector<Energy> _massa,_massb;

  /**
   * The function \f$\frac{\\hat{H}}{dq^2}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<double> _dh;

  /**
   * The function \f$\\hat{H}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<Energy2> _hres;

  /**
   * The \f$H(0)\f$ parameter  for the GS form of the
   *  Breit-Wigner
   */
  vector<Energy2> _h0;
  //@}
};

}


#endif /* HERWIG_TwoPionRhoCurrent_H */
