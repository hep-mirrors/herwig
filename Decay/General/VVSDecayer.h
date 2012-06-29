// -*- C++ -*-
#ifndef THEPEG_VVSDecayer_H
#define THEPEG_VVSDecayer_H
//
// This is the declaration of the VVSDecayer class.
//

#include "GeneralTwoBodyDecayer.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VVSVertexPtr;

/** \ingroup Decay
 * The VVSDecayer class implements the decay of a vector to a
 * vector and a scalar in a general model. It holds an VVSVertex pointer
 * that must be typecast from the  VertexBase pointer helid in the
 * GeneralTwoBodyDecayer. It implents the virtual functions me2() and
 * partialWidth(). 
 *
 * @see \ref VVSDecayerInterfaces "The interfaces"
 * defined for VVSDecayer.
 */
class VVSDecayer: public GeneralTwoBodyDecayer {

public:

  /**
   * The default constructor.
   */
  VVSDecayer() {}

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Return the matrix element squared for a given mode and phase-space channel
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
   * Initialize this object after the setup phase before saving an
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
  static ClassDescription<VVSDecayer> initVVSDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VVSDecayer & operator=(const VVSDecayer &);

private:

  /**
   *  Abstract pointer to AbstractVVSVertex
   */
  AbstractVVSVertexPtr _abstractVertex;

  /**
   * Pointer to the perturbative vertex
   */
  VVSVertexPtr _perturbativeVertex;

  /**
   *  Spin density matrix
   */
  mutable RhoDMatrix _rho;

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
 *  base classes of VVSDecayer. */
template <>
struct BaseClassTrait<VVSDecayer,1> {
  /** Typedef of the first base class of VVSDecayer. */
  typedef Herwig::GeneralTwoBodyDecayer NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VVSDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<VVSDecayer>
  : public ClassTraitsBase<VVSDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VVSDecayer"; }
};

/** @endcond */

}


#endif /* THEPEG_VVSDecayer_H */
