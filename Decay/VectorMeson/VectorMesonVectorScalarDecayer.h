// -*- C++ -*-
//
// VectorMesonVectorScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VectorMesonVectorScalarDecayer_H
#define HERWIG_VectorMesonVectorScalarDecayer_H
//
// This is the declaration of the VectorMesonVectorScalarDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace Herwig;

/** \ingroup Decay
 *
 *  This class is designed for the decay of a vector meson to a spin-1
 *  particle, either a vector meson or a photon, and a scalar meson.
 *  The current for the decay is
 *
 *  \f[\mathcal{M}=g\epsilon_mu\left[ p_V \cdot p_0 \epsilon_V^\mu  
 *                     -p_V^\mu \epsilon_V \cdot p_0\right],\f]
 *  where \f$\epsilon\f$ is the polarization vector of the decaying vector
 *  meson, \f$\epsilon_V\f$ is the polarization vector of the outgoing vector meson,
 *  \f$p_0\f$ is the momentum of the decaying particle and \f$p_V\f$ is the momentum
 *  of the outgoing vector meson.
 *
 * @see DecayIntegrator 
 * @see \ref VectorMesonVectorScalarDecayerInterfaces "The interfaces"
 * defined for VectorMesonVectorScalarDecayer.
 * 
 *  \author Peter Richardson
 * 
 */
class VectorMesonVectorScalarDecayer: public DecayIntegrator {

public:

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

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 4.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  VectorMesonVectorScalarDecayer & operator=(const VectorMesonVectorScalarDecayer &) = delete;

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
   * PDG codes for the incoming particles
   */
  vector<int> incoming_;

  /**
   * PDG codes for the outgoing particles (vector,scalar)
   */
  vector<pair<int,int> > outgoing_;

  /**
   * maximum weight for a decay
   */
  vector<double> maxWeight_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> vectors_[2];

};

}


#endif /* HERWIG_VectorMesonVectorScalarDecayer_H */
