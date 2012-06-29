// -*- C++ -*-
//
// UEDW0W1W1Vertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDW0W1W1Vertex_H
#define HERWIG_UEDW0W1W1Vertex_H
//
// This is the declaration of the UEDW0W1W1Vertex class.
//

#include "ThePEG/Helicity/Vertex/Vector/VVVVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::Direction;

/**
 * The coupling for the \f$ W\,W^{(1)}Z^{(1)}\f$ and  \f$ W\,W^{(1)}\gamma^{(1)}\f$ 
 * vertex . The class inherits from VVVVertex and implements 
 * the setCoupling member.
 *
 * @see \ref UEDW0W1W1VertexInterfaces "The interfaces"
 * defined for UEDW0W1W1Vertex.
 */
class UEDW0W1W1Vertex: public VVVVertex {

public:

  /**
   * The default constructor.
   */
  UEDW0W1W1Vertex();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<UEDW0W1W1Vertex> initUEDW0W1W1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDW0W1W1Vertex & operator=(const UEDW0W1W1Vertex &);

private:

  /**
   * The value of \f$\sin\theta_W\f$.
   */
  double theSinW;

  /**
   * The value of \f$\cos\theta_W\f$.
   */
  double theCosW;

  /**
   * The value of \f$\sin\theta_1\f$.
   */
  double theSinThetaOne;

  /**
   * The value of \f$\cos\theta_1\f$.
   */
  double theCosThetaOne;

  /**
   * The scale at which the coupling was last evaluated.
   */
  Energy2 theq2last;

  /**
   * The value of \f$ \sqrt{4\pi\alpha}\f$ when it was last evaluated.
   */
  Complex theElast;

  /**
   * The value of the coupling when it was last evaluated.
   */
  Complex theCouplast;

  /**
   * The PDG code of the last SM particle in the vertex 
   */
  long theSMlast;

  /**
   * The PDG code of the last KK-particle in the vertex 
   */
  long theKKlast;
};
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDW0W1W1Vertex. */
template <>
struct BaseClassTrait<Herwig::UEDW0W1W1Vertex,1> {
  /** Typedef of the first base class of UEDW0W1W1Vertex. */
  typedef ThePEG::Helicity::VVVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDW0W1W1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEDW0W1W1Vertex>
  : public ClassTraitsBase<Herwig::UEDW0W1W1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEDW0W1W1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDW0W1W1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDW0W1W1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEDW0W1W1Vertex_H */
