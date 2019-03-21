// -*- C++ -*-
//
// UEDZ0A1h1Vertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDZ0A1h1Vertex_H
#define HERWIG_UEDZ0A1h1Vertex_H
//
// This is the declaration of the UEDZ0A1h1Vertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * This class implements the coupling for the \f$A_{(1)}h_{(1)}Z^\mu_{{0}}\f$
 * vertex. It inherits from VSSVertex and implements the setCoupling member.
 *
 * @see \ref UEDZ0A1h1VertexInterfaces "The interfaces"
 * defined for UEDZ0A1h1Vertex.
 */
class UEDZ0A1h1Vertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  UEDZ0A1h1Vertex();

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
  static ClassDescription<UEDZ0A1h1Vertex> initUEDZ0A1h1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDZ0A1h1Vertex & operator=(const UEDZ0A1h1Vertex &) = delete;

private:
  
  /**
   * The value of \f$ \sin 2\theta_W \f$. 
   */
  double theSin2ThetaW;

  /**
   * The value of \f$m_kk/\sqrt{m_kk^2 + m_z^2}\f$
   */
  double theKappa;
    
  /**
   * The scale at which the coupling was last evaluated. 
   */
  Energy2 theq2Last;
  
  /**
   * The value of the coupling when it was last evaluated.
   */
  Complex theCoupLast;  
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDZ0A1h1Vertex. */
template <>
struct BaseClassTrait<Herwig::UEDZ0A1h1Vertex,1> {
  /** Typedef of the first base class of UEDZ0A1h1Vertex. */
  typedef ThePEG::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDZ0A1h1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEDZ0A1h1Vertex>
  : public ClassTraitsBase<Herwig::UEDZ0A1h1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEDZ0A1h1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDZ0A1h1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDZ0A1h1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEDZ0A1h1Vertex_H */
