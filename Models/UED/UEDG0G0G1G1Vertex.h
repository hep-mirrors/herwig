// -*- C++ -*-
//
// UEDG0G0G1G1Vertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDG0G0G1G1Vertex_H
#define HERWIG_UEDG0G0G1G1Vertex_H
//
// This is the declaration of the UEDG0G0G1G1Vertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVVertex.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the coupling of two Standard 
 * Model gluons to a pair of level 1 KK gluons.
 *
 * @see \ref UEDG0G0G1G1VertexInterfaces "The interfaces"
 * defined for UEDG0G0G1G1Vertex.
 */
class UEDG0G0G1G1Vertex: public Helicity::VVVVVertex {

public:

  /**
   * The default constructor.
   */
  UEDG0G0G1G1Vertex();

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
   *@param part4 The fourth interacting particle 
   */
  virtual void setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3, tcPDPtr part4);

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
  UEDG0G0G1G1Vertex & operator=(const UEDG0G0G1G1Vertex &);

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

#endif /* HERWIG_UEDG0G0G1G1Vertex_H */
