// -*- C++ -*-
//
// SMHiggsFermionsDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsFermionsDecayer_H
#define HERWIG_SMHiggsFermionsDecayer_H
//
// This is the declaration of the SMHiggsFermionsDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The SMHiggsFermionsDecayer class is designed to decay the Standard Model Higgs
 * to the Standard Model fermions.
 *
 * @see DecayIntegrator
 */
class SMHiggsFermionsDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  SMHiggsFermionsDecayer();
  
  /**
   * Which of the possible decays is required
   */
  virtual int modeNumber(bool & , tcPDPtr , const tPDVector & ) const {return -1;}

  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
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
  static ClassDescription<SMHiggsFermionsDecayer> initSMHiggsFermionsDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsFermionsDecayer & operator=(const SMHiggsFermionsDecayer &);

private:

  /**
   * Pointer to the Higgs vertex
   */
  AbstractFFSVertexPtr _hvertex;

  /**
   * maximum weights for the different decay modes
   */
  vector<double> _maxwgt;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   * Scalar wavefunction
   */
  mutable ScalarWaveFunction _swave;

  /**
   *  Spinor wavefunction
   */
  mutable vector<SpinorWaveFunction> _wave;

  /**
   *  Barred spinor wavefunction
   */
  mutable vector<SpinorBarWaveFunction> _wavebar;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsFermionsDecayer. */
template <>
struct BaseClassTrait<Herwig::SMHiggsFermionsDecayer,1> {
  /** Typedef of the first base class of SMHiggsFermionsDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsFermionsDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsFermionsDecayer>
  : public ClassTraitsBase<Herwig::SMHiggsFermionsDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsFermionsDecayer"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SMHiggsFermionsDecayer class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwPerturbativeHiggsDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHiggsFermionsDecayer_H */
