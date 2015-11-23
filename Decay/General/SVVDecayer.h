// -*- C++ -*-
//
// SVVDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SVVDecayer_H
#define HERWIG_SVVDecayer_H
//
// This is the declaration of the SVVDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.fh"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.fh"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VVSVertexPtr; 

/** \ingroup Decay
 * This SVVDecayer class implements the decay of a scalar to 
 * 2 vector bosons using either the tree level VVSVertex or the loop vertex.
 * It inherits from 
 * GeneralTwoBodyDecayer and implements the virtual member functions me2() 
 * and partialWidth(). It also stores a pointer to the VVSVertex.
 *
 * @see GeneralTwoBodyDecayer 
 * 
 */
class SVVDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  SVVDecayer() {}

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
  static ClassDescription<SVVDecayer> initSVVDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SVVDecayer & operator=(const SVVDecayer &);

private:
  
  /**
   *  Abstract pointer to general VVS vertex
   */
  AbstractVVSVertexPtr _abstractVertex;

  /**
   * Pointer to the perturbative form
   */
  VVSVertexPtr _perturbativeVertex; 

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Scalar wavefunction
   */
  mutable Helicity::ScalarWaveFunction _swave;

  /**
   *  Vector wavefunctions
   */
  mutable vector<Helicity::VectorWaveFunction> _vectors[2];
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SVVDecayer. */
template <>
struct BaseClassTrait<Herwig::SVVDecayer,1> {
  /** Typedef of the first base class of SVVDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SVVDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SVVDecayer>
  : public ClassTraitsBase<Herwig::SVVDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SVVDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_SVVDecayer_H */
