// -*- C++ -*-
//
// UEDF1F1G0Vertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDF1F1G0Vertex_H
#define HERWIG_UEDF1F1G0Vertex_H
//
// This is the declaration of the UEDF1F1G0Vertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This is the implementation of the \f$\bar{f}^{1}f^{(1)}g\f$ vertex
 *
 * @see \ref UEDF1F1G0VertexInterfaces "The interfaces"
 * defined for UEDF1F1G0Vertex.
 */
class UEDF1F1G0Vertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  UEDF1F1G0Vertex();

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
  virtual void setCoupling(Energy2 q2, tcPDPtr part1, tcPDPtr part2,
			   tcPDPtr part3);
  
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<UEDF1F1G0Vertex> initUEDF1F1G0Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDF1F1G0Vertex & operator=(const UEDF1F1G0Vertex &);

private:

  /**
   * The energy at which the coupling was last evaluated 
   */
  Energy2 theq2Last;

  /**
   * The value of the coupling when it was last evaluated 
   */
  Complex theCoupLast;
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDF1F1G0Vertex. */
template <>
struct BaseClassTrait<Herwig::UEDF1F1G0Vertex,1> {
  /** Typedef of the first base class of UEDF1F1G0Vertex. */
  typedef ThePEG::Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDF1F1G0Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEDF1F1G0Vertex>
  : public ClassTraitsBase<Herwig::UEDF1F1G0Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEDF1F1G0Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDF1F1G0Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDF1F1G0Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEDF1F1G0Vertex_H */
