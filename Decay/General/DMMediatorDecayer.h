// -*- C++ -*-
//
// DMMediatorDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DMMediatorDecayer_H
#define HERWIG_DMMediatorDecayer_H
//
// This is the declaration of the DMMediatorDecayer class.
//

#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "DMMediatorDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::VectorWaveFunction;

/**
 * Here is the documentation of the DMMediatorDecayer class.
 *
 * @see \ref DMMediatorDecayerInterfaces "The interfaces"
 * defined for DMMediatorDecayer.
 */
class DMMediatorDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  DMMediatorDecayer() : wgtmax_(0.) {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent,const tPDVector & children) const;
  //@}
  
  /** @name Virtual functions required by the Decayer class. */
  //@{

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
   * Function to return partial Width
   * @param inpart Pointer to incoming particle data object
   * @param out The outgoing particles from the current
   */
  virtual Energy partialWidth(tPDPtr inpart, vector<tPDPtr> out);
  //@}

  /**
   *  set up the decay
   */
  void setDecayInfo(PDPtr in, const vector<tPDPtr> & outCurrent,
		    WeakCurrentPtr current,
		    double normin, vector<double> sm);

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  The number of the mode
   * @param cc Whether of not this is the charge conjugate of the defined mode
   * @param id The PDG codes of the particles
   */
  int modeNumber(bool & cc, vector<long> id) const;

  /**
   *  Access to the map between the number of the mode and the modes in
   *  the current
   */
  unsigned int mode() const { return mode_; }

  /**
   *  Access to the weak current
   */
  WeakCurrentPtr weakCurrent() const { return current_; }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DMMediatorDecayer & operator=(const DMMediatorDecayer &) = delete;

private:
  
  /**
   * Incoming particle
   **/
  PDPtr inpart_;

  /**
   *  Outgoing particles from the current
   */
  vector<tPDPtr> currentOut_;
  
  /**
   * DM coupling to the dark mediator
   */
  double cDMmed_;

  /**                                                                                                                                                       
   * SM couplings to the dark mediator
   */
  vector<double> cSMmed_;

  /**
   * Pointer to the current
   */
  WeakCurrentPtr current_;

  /**
   * mapping of the modes to the currents
   */
  unsigned int mode_;

  /**
   * location of the weights
   */
  int wgtloc_;

  /**
   * the maximum weight
   */
  double wgtmax_;

  /**
   *  The weights for the different channels
   */
  vector<double> weights_;

  /**
   *  Spin density matrix 
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors for the decaying particle
   */
  mutable vector<VectorWaveFunction> vectors_;
};

}

#endif /* HERWIG_DMMediatorDecayer_H */
