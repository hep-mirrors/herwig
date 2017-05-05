// -*- C++ -*-
//
// TFFDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TFFDecayer_H
#define HERWIG_TFFDecayer_H
//
// This is the declaration of the TFFDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/Tensor/FFTVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/FFVTVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::FFTVertexPtr;

/** \ingroup Decay
 * The TFFDecayer class implements the decay of a tensor
 * to 2 fermions in a general model. It holds an FFTVertex pointer
 * that must be typecast from the VertexBase pointer held in 
 * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
 * partialWidth(). 
 *
 * @see GeneralTwoBodyDecayer
 */
class TFFDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  TFFDecayer() {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for.
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the matrix element
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

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay,MEOption meopt);
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
  static ClassDescription<TFFDecayer> initTFFDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TFFDecayer & operator=(const TFFDecayer &);

private:

  /**
   *  Abstract pointer to AbstractFFTVertex
   */
  AbstractFFTVertexPtr _abstractVertex;

  /**
   * Pointer to the perturbative vertex
   */
  FFTVertexPtr _perturbativeVertex;

  /**
   *  Abstract pointer to AbstractFFVVertex for QCD radiation from outgoing (anti)fermion
   */
  AbstractFFVVertexPtr _abstractOutgoingVertex1;

  /**
   *  Abstract pointer to AbstractFFVVertex for QCD radiation from outgoing (anti)fermion
   */
  AbstractFFVVertexPtr _abstractOutgoingVertex2;

  /**
   *  Abstract pointer to AbstractFFVTVertex for QCD radiation from 4 point vertex
   */
  AbstractFFVTVertexPtr _abstractFourPointVertex;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization tensors for the decaying particle
   */
  mutable vector<TensorWaveFunction> _tensors;

  /**
   *  Spinors for the decay products
   */
  mutable vector<SpinorWaveFunction> _wave;

  /**
   *  Barred spinors for the decay products
   */
  mutable vector<SpinorBarWaveFunction> _wavebar;

  /**
   *  Spin density matrix for 3 body decay
   */
  mutable RhoDMatrix _rho3;

  /**
   *  Tensor wavefunction for 3 body decay
   */
  mutable vector<TensorWaveFunction> _tensors3;

  /**
   *  Spinor wavefunction for 3 body decay
   */
  mutable vector<SpinorWaveFunction> _wave3;

  /**
   *  Barred spinor wavefunction for 3 body decay
   */
  mutable vector<SpinorBarWaveFunction> _wavebar3;

    /**
   *  Vector wavefunction for 3 body decay
   */
  mutable vector<VectorWaveFunction> _gluon;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TFFDecayer. */
template <>
struct BaseClassTrait<Herwig::TFFDecayer,1> {
  /** Typedef of the first base class of TFFDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TFFDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TFFDecayer>
  : public ClassTraitsBase<Herwig::TFFDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TFFDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_TFFDecayer_H */
