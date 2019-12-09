// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LeptoquarkModelSLQSLQGVertex_H
#define HERWIG_LeptoquarkModelSLQSLQGVertex_H
//
// This is the declaration of the LeptoquarkModelSLQSLQGVertex class.

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig/Models/Leptoquarks/LeptoquarkModel.h"

namespace Herwig {
using namespace ThePEG;
    
/** \ingroup Helicity
 * 
 *  The LeptoquarkModelSLQSLQGVertex class is the implementation of the gluon
 *  coupling to pairs of scalar leptoquarks in the LeptoquarkModel. It inherits from the VSSVertex 
 *  and implements the setCoupling member
 *
 *  @see VSSVertex
 *  @see VertexBase
 */
class LeptoquarkModelSLQSLQGVertex: public VSSVertex {
  
public:
  
  /**
   * Default constructor.
   */
  LeptoquarkModelSLQSLQGVertex();
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
  /**
   * Calculate the couplings. 
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3);

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
  
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
    
private:
  
  /**
   * Private and non-existent assignment operator.
   */
  LeptoquarkModelSLQSLQGVertex & operator=(const LeptoquarkModelSLQSLQGVertex &) = delete;
  
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

#endif /* HERWIG_LeptoquarkModelSLQSLQGVertex_H */
