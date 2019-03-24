// -*- C++ -*-
#ifndef HERWIG_SextetPSSVertex_H
#define HERWIG_SextetPSSVertex_H
//
// This is the declaration of the SextetPSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SextetPSSVertex class.
 *
 * @see \ref SextetPSSVertexInterfaces "The interfaces"
 * defined for SextetPSSVertex.
 */
class SextetPSSVertex: public Helicity::VSSVertex {

public:

  /**
   * The default constructor.
   */
  SextetPSSVertex() {
    colourStructure(ColourStructure::DELTA);
  }

 /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);

public:

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SextetPSSVertex & operator=(const SextetPSSVertex &) = delete;

private:
  
  /**
   * Store the value of the coupling when last evaluated
   */
  Complex coupLast_;

  /**
   * Store the scale at which coupling was last evaluated
   */
  Energy2 q2Last_;

};

}

#endif /* HERWIG_SextetPSSVertex_H */
