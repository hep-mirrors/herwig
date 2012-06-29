// -*- C++ -*-
//
// ScalarVectorVectorDecayer.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarVectorVectorDecayer_H
#define HERWIG_ScalarVectorVectorDecayer_H
//
// This is the declaration of the ScalarVectorVectorDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ScalarVectorVectorDecayer class is designed to perform the decay of 
 * a scalar meson to two spin-1 particles. The matrix element is taken
 * to have the form
 *  \f[\mathcal{M}=g\left[ p_1 \cdot p_2 \epsilon_1 \cdot \epsilon_2
 *                        -p_1 \cdot \epsilon_2 p_2 \cdot\epsilon_1\right],\f]
 *  where \f$\epsilon_{1,2}\f$ are the polarzation
 *  vectors of the outgoing vector particles.
 *
 * The incoming scalar meson, the outgoing vectors and the coupling can
 * be specified using the relevant interfaces.
 *
 * @see DecayIntegrator
 */
class ScalarVectorVectorDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  ScalarVectorVectorDecayer();
  
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
  double me2( const int ichan,const Particle & part,
	     const ParticleVector & decay, MEOption meopt) const;

  /**
   * Specify the \f$1\to2\f$ matrix element to be used in the running width calculation.
   * @param dm The DecayMode
   * @param mecode The code for the matrix element as described
   *               in the GenericWidthGenerator class, in this case 3.
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ScalarVectorVectorDecayer> initScalarVectorVectorDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ScalarVectorVectorDecayer & operator=(const ScalarVectorVectorDecayer &);

private:

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * the PDG code for the first outgoing particle
   */
  vector<int> _outgoing1;

  /**
   * the PDG code for the second outgoing particle
   */
  vector<int> _outgoing2;

  /**
   * the coupling for the decay, \f$g\f$.
   */
  vector<InvEnergy> _coupling;

  /**
   * the maximum weight for the decay
   */
  vector<double> _maxweight;

  /**
   *  initial number of modes
   */
  unsigned int _initsize;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the ecay products
   */
  mutable vector<Helicity::LorentzPolarizationVector> _vectors[2];
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ScalarVectorVectorDecayer. */
template <>
struct BaseClassTrait<Herwig::ScalarVectorVectorDecayer,1> {
  /** Typedef of the first base class of ScalarVectorVectorDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ScalarVectorVectorDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ScalarVectorVectorDecayer>
  : public ClassTraitsBase<Herwig::ScalarVectorVectorDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ScalarVectorVectorDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the ScalarVectorVectorDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwSMDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ScalarVectorVectorDecayer_H */
