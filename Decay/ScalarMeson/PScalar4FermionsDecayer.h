// -*- C++ -*-
//
// PScalar4FermionsDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PScalar4FermionsDecayer_H
#define HERWIG_PScalar4FermionsDecayer_H
// This is the declaration of the PScalar4FermionsDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzSpinorBar.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>PScalar4FermionsDecayer</code> class is designed for the
 *  decay of the neutral pion to four leptons, in this case electrons.
 *  The propagator for the off-shell boson is taken to be \f$\frac1{m^2}\f$. There
 *  is also the option of including a vector meson dominance  form-factor.
 *
 *  In this case the matrix element is
 *  \f[\mathcal{M} = \frac{g}{m^2_{f_1\bar{f_1}}m^2_{f_2\bar{f_2}}}
 *                   \epsilon^{\mu\nu\alpha\beta}    
 *                   \bar{u}(p_{f_1})\gamma_\mu    v(p_{\bar{f_1}})p_{f_1\bar{f_1}\nu}
 *                   \bar{u}(p_{f_2})\gamma_\alpha v(p_{\bar{f_2}}) p_{f_2\bar{f_2}\beta}
 *  \f]
 *  It includes the option of a vector meson dominance (VMD) type form factor  
 *  \f$\frac{-M^2+i\Gamma M}{(m^2_{f\bar{f}}-M^2+i\Gamma M)}\f$. In the case of identical
 *  fermions in also includes the exchange diagram.
 *  
 * @see DecayIntegrator
 * @see PScalarVectorVectorDecayer
 * @see PScalarVectorFermionsDecayer
 * 
 *  \author Peter Richardson
 *
 */
class PScalar4FermionsDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  PScalar4FermionsDecayer();
  
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
  PScalar4FermionsDecayer & operator=(const PScalar4FermionsDecayer &) = delete;

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
   * the PDG codes for the incoming particle
   */
  vector<int> incoming_;

  /**
   * the PDG codes for the first outgoing fermion
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
   * PDG code of the particle for the VMD factor
   */
  vector<int> VMDid_;

  /**
   * Mass of the particle for the VMD factor
   */
  vector<Energy> VMDmass_;

  /**
   * Width of the particle for the VMD factor
   */
  vector<Energy> VMDwidth_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinor   <SqrtEnergy> > wave_[2];

  /**
   *  Barred spinors for the decay products
   */
  mutable vector<Helicity::LorentzSpinorBar<SqrtEnergy> > wavebar_[2];
};

}


#endif /* HERWIG_PScalar4FermionsDecayer_H */
