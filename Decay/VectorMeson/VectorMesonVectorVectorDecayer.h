// -*- C++ -*-
//
// VectorMesonVectorVectorDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_VectorMesonVectorVectorDecayer_H
#define HERWIG_VectorMesonVectorVectorDecayer_H
//
// This is the declaration of the VectorMesonVectorVectorDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {
using namespace  ThePEG;

/** \ingroup Decay
 *
 *  This class is designed for the decay of a vector meson to two spin
 *  one particles, either other vector mesons or photons.
 *
 *  The current for the decay is taken to be
 *
 *  \f[\mathcal{M}= \frac{g}{M_0^2}
 *             ( p_{0\nu}\epsilon^\alpha-p_{0\alpha} \epsilon^\nu)\left[
 *             ((p_{1\nu} \epsilon_1^\beta- p_1^\beta \epsilon_{1\nu})
 *             (p_{2\alpha} \epsilon_{2\beta}- p_{2\beta} \epsilon_{2\alpha})
 *             -(\nu \leftrightarrow\alpha))\right],\f]
 * where \f$p_{0,1,2}\f$ are the momenta of the incoming and two outgoing
 * vector mesons and \f$\epsilon_{1,2}\f$ are the polarization vectors of
 * the outgoing mesons and \f$\epsilon\f$ is the polarization vector of the 
 * incoming meson. 
 *
 *  The incoming vector mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$V\to VV\f$ decays are specified in the default
 *  constructor.                                 
 *
 * @see DecayIntegrator.
 * @see \ref VectorMesonVectorVectorDecayerInterfaces "The interfaces"
 * defined for VectorMesonVectorVectorDecayer.
 */
class VectorMesonVectorVectorDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  VectorMesonVectorVectorDecayer();

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 5.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

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
  VectorMesonVectorVectorDecayer & operator=(const VectorMesonVectorVectorDecayer &) = delete;

private:

  /**
   * the coupling for the decay
   */
  vector<double> _coupling;

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * the PDG codes for the 1st outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG codes for the 2nd outgoing particle
   */
  vector<int> _outgoing2;

  /**
   * maximum weight for a decay
   */
  vector<double> _maxweight;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors[3];


};

}


#endif /* HERWIG_VectorMesonVectorVectorDecayer_H */
