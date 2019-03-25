// -*- C++ -*-
#ifndef HERWIG_SextetGGSSVertex_H
#define HERWIG_SextetGGSSVertex_H
//
// This is the declaration of the SextetGGSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VVSSVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SextetGGSSVertex class.
 *
 * @see \ref SextetGGSSVertexInterfaces "The interfaces"
 * defined for SextetGGSSVertex.
 */
class SextetGGSSVertex: public Helicity::VVSSVertex {

public:

  /**
   * The default constructor.
   */
  SextetGGSSVertex();

  /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle. 
  */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
 			   tcPDPtr part2,tcPDPtr part3,tcPDPtr part4);

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
  SextetGGSSVertex & operator=(const SextetGGSSVertex &) = delete;

private:

  /**
   * The energy at which the coupling was last evaluated
   */
  Energy2 q2last_;

  /**
   * The coupling when it was last evaluated
   */
  Complex couplast_;

};

}

#endif /* HERWIG_SextetGGSSVertex_H */
