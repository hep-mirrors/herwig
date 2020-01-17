// -*- C++ -*-
//
// UEDG1G1G0Vertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDG1G1G0Vertex_H
#define HERWIG_UEDG1G1G0Vertex_H
//
// This is the declaration of the UEDG1G1G0Vertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::Direction;

/**
 * This is the implementation of the Feynman rule for the coupling
 * of Standard Model gluon to a pair of (level1) KK excited gluons. 
 *
 * @see \ref UEDG1G1G0VertexInterfaces "The interfaces"
 * defined for UEDG1G1G0Vertex.
 */
class UEDG1G1G0Vertex: public VVVVertex {

public:

  /**
   * The default constructor.
   */
  UEDG1G1G0Vertex();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

  /** Calculate the coupling
   *@param q2 The scale at which to evaluate the coupling
   *@param part1 The first interacting particle 
   *@param part2 The second interacting particle 
   *@param part3 The third interacting particle 
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
  UEDG1G1G0Vertex & operator=(const UEDG1G1G0Vertex &) = delete;

private:
  
  /**
   * The scale at which the coupling was last evaluated
   */
  Energy2 theq2Last;
  
  /**
   * The value of the coupling when it was last evaluated
   */
  Complex theCoupLast;
  
};
}

#endif /* HERWIG_UEDG1G1G0Vertex_H */
