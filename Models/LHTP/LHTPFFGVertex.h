// -*- C++ -*-
#ifndef HERWIG_LHTPFFGVertex_H
#define HERWIG_LHTPFFGVertex_H
//
// This is the declaration of the LHTPFFGVertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The LHTPFFGVertex class implements the coupling of the 
 * gluon to the coloured fermions, the SM quarks, the extra top-like quark
 * and the T-parity odd quarks of the Little Higgs model with T-parity.
 */
class LHTPFFGVertex: public Helicity::FFVVertex {

public:

  /**
   * The default constructor.
   */
  LHTPFFGVertex();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  LHTPFFGVertex & operator=(const LHTPFFGVertex &) = delete;

private:

  /**
   * Storage of the couplings.
   */
  //@{
  /**
   *  The last value of the strong coupling calculated.
   */
  Complex coupLast_;

  /**
   *  The scale \f$q^2\f$ at which the coupling was last evaluated.
   */
  Energy2 q2Last_;
  //@}

};

}

#endif /* HERWIG_LHTPFFGVertex_H */
