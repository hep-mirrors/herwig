// -*- C++ -*-
//
// PScalarVectorFermionsDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_PScalarVectorFermionsDecayer_H
#define THEPEG_PScalarVectorFermionsDecayer_H
//
// This is the declaration of the PScalarVectorFermionsDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 * The <code>PScalarVectorFermionsDecayer</code> class is designed for the decay of a 
 * pseudoscalar meson to a spin-1 particle and a fermion-antifermion pair. In practice
 * these decays are of the form \f$\gamma\ell^+\ell^-\f$ and the propagator of
 * the off-shell boson is taken to be \f$\frac1{m^2_{f\bar{f}}}\f$.
 * There is also the option of including a vector meson dominance
 * form-factor.
 *
 *  In this case the matrix element is
 *  \f[\mathcal{M} = \frac{g}{m^2_{f\bar{f}}}
 *                   \epsilon^{\mu\nu\alpha\beta}p_{V\mu}\epsilon_{V\nu}
 *                   \bar{u}(p_f)\gamma_\alpha v(p_{\bar{f}}) p_{f\bar{f}\beta}
 *  \f]
 *  It includes the option of a vector meson dominance (VMD) type form factor  
 *  \f$\frac{-M^2+i\Gamma M}{(m^2_{f\bar{f}}-M^2+i\Gamma M)}\f$.
 *    
 *  The incoming pseudoscalar meson, the outgoing vector, the fermion and antifermion
 *  and the coupling can be specified using the relevant interfaces.
 *
 * @see DecayIntegrator
 * @see PScalarVectorVectorDecayer
 * @see PScalar4FermionsDecayer
 * 
 *  \author Peter Richardson
 *
 */
class PScalarVectorFermionsDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  PScalarVectorFermionsDecayer();
  
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;
  
  /**
   * The differential three body decay rate with one integral performed.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s  The invariant mass which still needs to be integrate over.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The differential rate \f$\frac{d\Gamma}{ds}\f$
   */
  virtual InvEnergy threeBodydGammads(const int imode, const Energy2 q2, 
				      const Energy2 s,
				      const Energy m1, const Energy m2, 
				      const Energy m3) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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

  /**
   * Initialize this object to the begining of the run phase.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Private and non-existent assignment operator.
   */
  PScalarVectorFermionsDecayer & operator=(const PScalarVectorFermionsDecayer &) = delete;

public:

  /**
   *   Set the parameters for a decay mode
   */
  string setUpDecayMode(string arg);

private:

  /**
   * coupling for a decay
   */
  vector<InvEnergy> coupling_;

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> incoming_;

  /**
   * the PDG codes for the outgoing vector
   */
  vector<pair<long,long> > outgoing_;

  /**
   * maximum weight for a decay
   */
  vector<double> maxweight_;

  /**
   * Include the VMD factor
   */
  vector<int> includeVMD_;

  /**
   * PDG code for thte particle to use in the VMD factor.
   */
  vector<int> VMDid_;

  /**
   * Mass to use in the VMD factor.
   */
  vector<Energy> VMDmass_;

  /**
   * Width to use in the VMD factor.
   */
  vector<Energy> VMDwidth_;
  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors for the decay product
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_;

  /**
   *  Spinors for the fermions
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_;

  /**
   *  Barred spinors for the fermions 
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_;
};

}


#endif /* THEPEG_PScalarVectorFermionsDecayer_H */
