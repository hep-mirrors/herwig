// -*- C++ -*-
//
// TensorMesonVectorVectorDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TensorMesonVectorVectorDecayer_H
#define HERWIG_TensorMesonVectorVectorDecayer_H
//
// This is the declaration of the TensorMesonVectorVectorDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzTensor.h"

namespace Herwig {
using namespace ThePEG; 

/** \ingroup Decay
 *
 *  The TensorMesonVectorVectorDecayer class is designed to simulate
 *  the decay of a tensor meson to two spin-1 particles. In practice, at least
 *  for the light tensor mesons, this is only the decay of a tensor meson to two
 *  photons. In principle for bottom and charm tensors this may be the decay to
 *  two vector mesons an d it is used to simulate three and four body 
 *  \f$a_2\f$, \f$f_2\f$ and \f$K_2\f$ decays. 
 *
 *  The form of the matrix element is based on the perturbative matrix element
 *  for the decay of a graviton to two vector bosons with the neglect of a mass term
 *
 * \f[ \mathcal{M} =  g\epsilon_{\mu\nu}\left[
 *   \left(\epsilon_{1\alpha} p_1^\mu - \epsilon_1^\mu p_{1\alpha}\right)
 *   \left(\epsilon_2^\alpha  p_2^\nu - \epsilon_2^\nu p_2^\alpha\right)
 *  +\left(\epsilon_{1\alpha} p_1^\nu - \epsilon_1^\nu p_{1\alpha}\right)
 *   \left(\epsilon_2^\alpha  p_2^\mu - \epsilon_2^\nu p_2^\alpha\right)
 * -\frac12g^{\mu\nu}
 *    \left(\epsilon_{1\alpha} p_{1\beta}- \epsilon_{1\beta} p_{1\alpha}\right)
 *    \left(\epsilon_2^\alpha  p_2^\beta - \epsilon_2^\beta p_2^\alpha\right)\right] \f]
 *  in such a way that it vanishes if the polarizations of the outgoing vectors are
 *  replaced with their momenta.
 *
 *  The incoming tensor mesons together with their decay products and the coupling 
 *  \f$g\f$ can be specified using the interfaces for the class. The maximum weights
 *  for the decays can be calculated using the Initialize interface of the
 *  DecayIntegrator class or specified using the interface.
 *
 *  The incoming and outgoing particles, couplings and maximum weights for
 *  many of the common \f$T\to VV\f$ decays are specified in the default
 *  constructor.
 *
 * @see DecayIntegrator
 * * @see \ref TensorMesonVectorVectorDecayerInterfaces "The interfaces"
 * defined for TensorMesonVectorVectorDecayer.
 *
 */
class TensorMesonVectorVectorDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  TensorMesonVectorVectorDecayer();

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
   *               in the GenericWidthGenerator class, in this case 9.
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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<TensorMesonVectorVectorDecayer> initTensorMesonVectorVectorDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  TensorMesonVectorVectorDecayer & operator=(const TensorMesonVectorVectorDecayer &) = delete;

private:

  /**
   * the PDG codes for the incoming particles
   */
  vector<int> _incoming;

  /**
   * the PDG codes for the first outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG codes for the second outgoing particle
   */
  vector<int> _outgoing2;

  /**
   * the coupling for the decay
   */
  vector<InvEnergy> _coupling;

  /**
   * the maximum weight for the decay
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
  mutable vector<Helicity::LorentzPolarizationVector > _vectors[2];

  /**
   *   Storage of the \f$\rho\f$ matrix
   */
  mutable RhoDMatrix _rho;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of Herwig::TensorMesonVectorVectorDecayer.
 */
template <>
struct BaseClassTrait<Herwig::TensorMesonVectorVectorDecayer,1> {
  /** Typedef of the base class of TensorMesonVectorVectorDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TensorMesonVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::TensorMesonVectorVectorDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::TensorMesonVectorVectorDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwTMDecay.so"; }
  
};
  
/** @endcond */
  
}

#endif /* HERWIG_TensorMesonVectorVectorDecayer_H */
