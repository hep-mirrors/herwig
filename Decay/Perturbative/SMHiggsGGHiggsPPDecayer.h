// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMHiggsGGHiggsPPDecayer_H
#define HERWIG_SMHiggsGGHiggsPPDecayer_H
//
// This is the declaration of the SMHiggsGGHiggsPPDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Models/StandardModel/SMHGGVertex.h"
#include "Herwig/Models/StandardModel/SMHPPVertex.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * Typedef for the \f$H\to gg\f$ vertex
 */
typedef Ptr<Herwig::SMHGGVertex>::pointer HGGPtr;

/**
 * Typedef for the \f$H\to \gamma\gamma\f$ vertex
 */
typedef Ptr<Herwig::SMHPPVertex>::pointer HPPPtr;
  
/**
 * The <code>SMHiggsGGHiggsPPDecayer</code> class performs the
 * of a Standard Model Higgs boson to either a pair
 * of photons or a pair of gluons.
 *
 * @see DecayIntegrator
 */ 
class SMHiggsGGHiggsPPDecayer: public DecayIntegrator {
  
public:

  /**
   * The default constructor.
   */
  SMHiggsGGHiggsPPDecayer() : _h0wgt(2,1.) {}
  
  /** @name Virtual functions required by the Decayer class. */
  //@{
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
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;
  
  /**
   * Which of the possible decays is required
   */
  virtual int modeNumber(bool &, tcPDPtr, const tPDVector & ) const {return -1;}
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,const tPDVector & children) const;
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
  static ClassDescription<SMHiggsGGHiggsPPDecayer> 
  initSMHiggsGGHiggsPPDecayer;
  
  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMHiggsGGHiggsPPDecayer & operator=(const SMHiggsGGHiggsPPDecayer &);
  
  /**
   * Pointer to h->gluon,gluon vertex
   */
  HGGPtr _hggvertex;
  
  /**
   * Pointer to h->gamma,gamma vertex
   */
  HPPPtr _hppvertex;
  
  /**
   * Maximum weight for integration
   */
  vector<double> _h0wgt;
  
  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Scalar wavefunction
   */
  mutable ScalarWaveFunction _swave;

  /**
   *  Vector wavefunctions
   */
  mutable vector<VectorWaveFunction> _vwave[2];
};
  
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMHiggsGGHiggsPPDecayer. */
template <>
struct BaseClassTrait<Herwig::SMHiggsGGHiggsPPDecayer,1> {
  /** Typedef of the first base class of SMHiggsGGHiggsPPDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMHiggsGGHiggsPPDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMHiggsGGHiggsPPDecayer>
  : public ClassTraitsBase<Herwig::SMHiggsGGHiggsPPDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMHiggsGGHiggsPPDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SMHiggsGGHiggsPPDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwPerturbativeHiggsDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMHiggsGGHiggsPPDecayer_H */
