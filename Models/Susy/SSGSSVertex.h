// -*- C++ -*-
//
// SSGSSVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SSGSSVertex_H
#define HERWIG_SSGSSVertex_H
//
// This is the declaration of the SSGSSVertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "SusyBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The SSGSSVertex implements the coupling of a gluon to 2 sfermions. 
 * It inherits from VSSVertex and implements the setCoupling method.
 * 
 * @see VSSVertex
 */
class SSGSSVertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  SSGSSVertex();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

 /**
   * Calculate the couplings.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
                           tcPDPtr part2,tcPDPtr part3);

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
  SSGSSVertex & operator=(const SSGSSVertex &);

private:
  
  /**
   * Store the value of the coupling when last evaluated
   */
  Complex _couplast;

  /**
   * Store the scale at which coupling was last evaluated
   */
  Energy2 _q2last;
};
}

#endif /* HERWIG_SSGSSVertex_H */
