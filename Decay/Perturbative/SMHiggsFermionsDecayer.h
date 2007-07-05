// -*- C++ -*-
#ifndef HERWIG_SMHiggsFermionsDecayer_H
#define HERWIG_SMHiggsFermionsDecayer_H
//
// This is the declaration of the SMHiggsFermionsDecayer class.
//

#include "Herwig++/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "SMHiggsFermionsDecayer.fh"

namespace Herwig {
using namespace ThePEG;
using Helicity::FFSVertexPtr;

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

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param dm The decay mode
   */
  virtual int modeNumber(bool &, const DecayMode &) const {return -1;}

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(bool vertex, const int ichan, const Particle & part,
		     const ParticleVector & decay) const;
  //@}

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
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
  virtual void doinit() throw(InitException);
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
  FFSVertexPtr _hvertex;

  /**
   * maximum weights for the different decay modes
   */
  vector<double> _maxwgt;

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
  static string className() { return "Herwig++::SMHiggsFermionsDecayer"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SMHiggsFermionsDecayer class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwPerturbativeHiggsDecay.so"; }
};

/** @endcond */

}

#include "SMHiggsFermionsDecayer.icc"

#endif /* HERWIG_SMHiggsFermionsDecayer_H */
