// -*- C++ -*-
//
// ScalarScalarScalarDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ScalarScalarScalarDecayer_H
#define HERWIG_ScalarScalarScalarDecayer_H
// This is the declaration of the ScalarScalarScalarDecayer class.

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>ScalarScalarScalarDecayer</code> class is designed for the decay
 *  of a scalar meson to two scalar mesons. This decay is obviously trival,
 *  however this decayer allows us to include off-shell effects for the outgoing
 *  particles and the correct calculation of the running width for the incoming
 *  particle.
 *
 *  The matrix element is taken to b
 *  \f[\mathcal{M} = g \f].
 *
 * @see DecayIntegrator
 * 
 * \author Peter Richardson
 *
 */
class ScalarScalarScalarDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  ScalarScalarScalarDecayer();
  
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
   *               in the GenericWidthGenerator class, in this case 6.
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
  static ClassDescription<ScalarScalarScalarDecayer> initScalarScalarScalarDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  ScalarScalarScalarDecayer & operator=(const ScalarScalarScalarDecayer &) = delete;

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
   * the PDG code for the second incoming particle
   */
  vector<int> _outgoing2;

  /**
   * the coupling for the decay
   */
  vector<Energy> _coupling;

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
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ScalarScalarScalarDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::ScalarScalarScalarDecayer,1> {
    /** Typedef of the base class of ScalarScalarScalarDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ScalarScalarScalarDecayer>
  : public ClassTraitsBase<Herwig::ScalarScalarScalarDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig::ScalarScalarScalarDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwSMDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_ScalarScalarScalarDecayer_H */
