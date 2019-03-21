// -*- C++ -*-
//
// SFFDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SFFDecayer_H
#define HERWIG_SFFDecayer_H
//
// This is the declaration of the SFFDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::FFSVertexPtr;

  /** \ingroup Decay
   * The SFFDecayer class implements the decay of a scalar to 2
   * fermions in a general model. It holds an FFSVertex pointer that 
   * must be typecast from the VertexBase pointer held in 
   * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
   * partialWidth(). 
   *
   * @see GeneralTwoBodyDecayer
   */
class SFFDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  SFFDecayer() {}

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

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay,MEOption meopt);

  /**
   * Indentify outgoing vertices for the fermion and antifermion
   */
  void identifyVertices(const int iferm, const int ianti,
			const Particle & inpart, const ParticleVector & decay,
			AbstractFFVVertexPtr & abstractOutgoingVertexF, 
			AbstractFFVVertexPtr & abstractOutgoingVertexA);
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
  static ClassDescription<SFFDecayer> initSFFDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SFFDecayer & operator=(const SFFDecayer &) = delete;

private:

  /**
   *  Abstract pointer to AbstractFFSVertex
   */
  AbstractFFSVertexPtr _abstractVertex;

  /**
   * Pointer to the perturbative vertex
   */
  FFSVertexPtr _perturbativeVertex;

  /**
   *  Abstract pointer to AbstractVSSVertex for QCD radiation from incoming scalar
   */
  AbstractVSSVertexPtr _abstractIncomingVertex;

  /**
   *  Abstract pointer to AbstractFFVVertex for QCD radiation from outgoing (anti)fermion
   */
  AbstractFFVVertexPtr _abstractOutgoingVertex1;

  /**
   *  Abstract pointer to AbstractFFVVertex for QCD radiation from outgoing (anti)fermion
   */
  AbstractFFVVertexPtr _abstractOutgoingVertex2;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

  /**
   *  Scalar wavefunction
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

 /**
   *  Spin density matrix for 3 body decay
   */
  mutable RhoDMatrix _rho3;

  /**
   *  Scalar wavefunction for 3 body decay
   */
  mutable ScalarWaveFunction _swave3;

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
 *  base classes of SFFDecayer. */
template <>
struct BaseClassTrait<Herwig::SFFDecayer,1> {
  /** Typedef of the first base class of SFFDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SFFDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SFFDecayer>
  : public ClassTraitsBase<Herwig::SFFDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SFFDecayer"; }
};

/** @endcond */

}


#endif /* HERWIG_SFFDecayer_H */
