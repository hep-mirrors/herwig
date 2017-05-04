// -*- C++ -*-
//
// TVVDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TVVDecayer_H
#define HERWIG_TVVDecayer_H
//
// This is the declaration of the TVVDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVTVertex.h"
#include "ThePEG/Helicity/Vertex/Tensor/VVVTVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VVTVertexPtr;

/** \ingroup Decay
 * The TVVDecayer class implements the decay of a tensor
 * to 2 vector bosons in a general model. It holds a VVTVertex pointer
 * that must be typecast from the VertexBase pointer held in 
 * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
 * partialWidth(). 
 *
 * @see GeneralTwoBodyDecayer
 */
class TVVDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  TVVDecayer() {}

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
  static ClassDescription<TVVDecayer> initTVVDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TVVDecayer & operator=(const TVVDecayer &);

private:

  /**
   *  Abstract pointer to AbstractVVTVertex
   */
  AbstractVVTVertexPtr _abstractVertex;

  /**
   * Pointer to the perturbative vertex
   */
  VVTVertexPtr _perturbativeVertex;

  /**
   *  Abstract pointer to AbstractVVVVertex for QCD radiation from outgoing vector
   */
  AbstractVVVVertexPtr _abstractOutgoingVertex1;

  /**
   *  Abstract pointer to AbstractVVVVertex for QCD radiation from outgoing vector
   */
  AbstractVVVVertexPtr _abstractOutgoingVertex2;

  /**
   *  Abstract pointer to AbstractVVVTVertex for QCD radiation from 4 point vertex
   */
  AbstractVVVTVertexPtr _abstractFourPointVertex;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization tensors of decaying particle
   */
  mutable vector<Helicity::TensorWaveFunction> _tensors;

  /**
   *  Polarization vectors of outgoing vector bosons
   */
  mutable vector<Helicity::VectorWaveFunction> _vectors[2];

  /**
   *  Spin density matrix for 3 body decay
   */
  mutable RhoDMatrix _rho3;

  /**
   *  Tensor wavefunction for 3 body decay
   */
  mutable vector<Helicity::TensorWaveFunction> _tensors3;

  /**
   *  Polarization vectors of outgoing vector bosons
   */
  mutable vector<Helicity::VectorWaveFunction> _vectors3[2];

    /**
   *  Vector wavefunction for 3 body decay
   */
  mutable vector<Helicity::VectorWaveFunction> _gluon;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TVVDecayer. */
template <>
struct BaseClassTrait<Herwig::TVVDecayer,1> {
  /** Typedef of the first base class of TVVDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TVVDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TVVDecayer>
  : public ClassTraitsBase<Herwig::TVVDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TVVDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_TVVDecayer_H */
