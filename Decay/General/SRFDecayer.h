// -*- C++ -*-
//
// SRFDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SRFDecayer_H
#define HERWIG_SRFDecayer_H
//
// This is the declaration of the SRFDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/Scalar/RFSVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::RFSVertexPtr;

/** \ingroup Decay
 * The SRFDecayer class implements the decay of a scalar to spin-3/2
 * and spin-1/2 fermion in a general model. It holds an RFSVertex pointer that 
 * must be typecast from the VertexBase pointer held in 
 * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
 * partialWidth(). 
 *
 * @see GeneralTwoBodyDecayer
 */
class SRFDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  SRFDecayer() {}

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
   * Function to return partial Width
   * @param inpart The decaying particle.
   * @param outa One of the decay products.
   * @param outb The other decay product.
   */
  virtual Energy partialWidth(PMPair inpart, PMPair outa, 
			      PMPair outb) const;
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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SRFDecayer> initSRFDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SRFDecayer & operator=(const SRFDecayer &);

private:

  /**
   *  Abstract pointer to AbstractFFSVertex
   */
  AbstractRFSVertexPtr abstractVertex_;

  /**
   * Pointer to the perturbative vertex
   */
  RFSVertexPtr perturbativeVertex_;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix rho_;

  /**
   *  Scalar wavefunction
   */
  mutable ScalarWaveFunction swave_;

  /**
   *  Spinor wavefunction
   */
  mutable vector<SpinorWaveFunction> wave_;

  /**
   *  Barred spinor wavefunction
   */
  mutable vector<SpinorBarWaveFunction> wavebar_;

  /**
   *  RS Spinor wavefunction
   */
  mutable vector<RSSpinorWaveFunction> RSwave_;

  /**
   *  Barred RS spinor wavefunction
   */
  mutable vector<RSSpinorBarWaveFunction> RSwavebar_;

  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SRFDecayer. */
template <>
struct BaseClassTrait<Herwig::SRFDecayer,1> {
  /** Typedef of the first base class of SRFDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SRFDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SRFDecayer>
  : public ClassTraitsBase<Herwig::SRFDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SRFDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_SRFDecayer_H */
