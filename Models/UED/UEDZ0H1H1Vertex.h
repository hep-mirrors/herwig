// -*- C++ -*-
//
// UEDZ0H1H1Vertex.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDZ0H1H1Vertex_H
#define HERWIG_UEDZ0H1H1Vertex_H
//
// This is the declaration of the UEDZ0H1H1Vertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/VSSVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The implementation of the \f$Z^\mu_{(0)}H_{(1)}^+H_{(1)}^-\f$ vertex. This
 * class inherits from VSSVertex and implements the setCoupling member function.
 * 
 * The vertex is taken to have the form:
 * \f[\frac{g}{1 + m_W^2R^2}\left(\frac{\cos 2\theta_W}{2\cos\theta_W} -
 * m_W^2 R^2\cos^2\theta_W \right)\left(p(H_{(1)}^-) - p(H_{(1)}^+)\right) \f]
 *
 * @see \ref UEDZ0H1H1VertexInterfaces "The interfaces"
 * defined for UEDZ0H1H1Vertex.
 */
class UEDZ0H1H1Vertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  UEDZ0H1H1Vertex();

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
  static ClassDescription<UEDZ0H1H1Vertex> initUEDZ0H1H1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDZ0H1H1Vertex & operator=(const UEDZ0H1H1Vertex &) = delete;

private:

  /**
   * The value of \f$\cos\theta_W\f$.
   */
  double theCosThetaW;
    
  /**
   * The value of \f$\cos 2\theta_W\f$.
   */
  double theCosTheta2W;

  /**
   * The mass-squared of the \f$W\f$ boson.
   */
  Energy2 theMw2;

  /**
   * The square of the compactification radius.  
   */
  InvEnergy2 theR2;  

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
 *  base classes of UEDZ0H1H1Vertex. */
template <>
struct BaseClassTrait<Herwig::UEDZ0H1H1Vertex,1> {
  /** Typedef of the first base class of UEDZ0H1H1Vertex. */
  typedef ThePEG::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDZ0H1H1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEDZ0H1H1Vertex>
  : public ClassTraitsBase<Herwig::UEDZ0H1H1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEDZ0H1H1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDZ0H1H1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDZ0H1H1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEDZ0H1H1Vertex_H */
