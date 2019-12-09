// -*- C++ -*-
//
// TensorMesonVectorPScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TensorMesonVectorPScalarDecayer_H
#define HERWIG_TensorMesonVectorPScalarDecayer_H
//
// This is the declaration of the TensorMesonVectorPScalarDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>TensorMesonVectorPScalarDecayer</code> class handles the decay of
 *  a tensor meson to a vector and a pseudoscalar. Examples of this decay
 *  are \f$a_2\to\rho\pi\f$ or \f$a_2 \to\pi \gamma\f$. The matrix element is
 *  given by
 *  \f[\mathcal{M}=\epsilon_T^{\mu\nu}p_{P,\mu} \epsilon_{\nu\alpha\beta\gamma}
 *     p_V^\alpha \epsilon_V^\beta p_P^\gamma,\f]
 * where \f$p_P\f$ is the momentum of the pseudoscalar meson, \f$p_V\f$ is
 * the momentum of the vector, \f$\epsilon_T\f$ is the polarization tensor of the 
 * decaying meson and \f$\epsilon_V\f$ is the polarization vector of the outgoing
 * vector meson.
 *
 *  The incoming tensor mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$T\to PV\f$ decays are specified in the default
 *  constructor.
 *
 * @see DecayIntegrator
 * @see \ref TensorMesonVectorPScalarDecayerInterfaces "The interfaces"
 * defined for TensorMesonVectorPScalarDecayer.
 */
class TensorMesonVectorPScalarDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  TensorMesonVectorPScalarDecayer();

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
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 8.
   * @param coupling The coupling for the matrix element.
   * @return True or False if this mode can be handled.
   */
  bool twoBodyMEcode(const DecayMode & dm, int & mecode, double & coupling) const;

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
  TensorMesonVectorPScalarDecayer & operator=(const TensorMesonVectorPScalarDecayer &) = delete;

private:

  /**
   * PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * PDG codes for the outgoing vector
   */
  vector<int> _outgoingV;

  /**
   * PDG codes for the outgoing pseudoscalar
   */
  vector<int> _outgoingP;

  /**
   * coupling for the decay
   */
  vector<InvEnergy2> _coupling;

  /**
   * max weight ofr the decay
   */
  vector<double> _maxweight;

  /**
   *  Initial size of the vectors
   */
  unsigned int _initsize;

  /**
   *  Storage of polarization tensors to try and increase
   *  speed
   */
  mutable vector<Helicity::LorentzTensor<double> > _tensors;

  /**
   *  Storage of the polarization vectors
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors;

  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix _rho;
};

}


#endif /* HERWIG_TensorMesonVectorPScalarDecayer_H */
