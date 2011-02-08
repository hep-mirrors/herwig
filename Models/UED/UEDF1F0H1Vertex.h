// -*- C++ -*-
//
// UEDF1F0H1Vertex.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_UEDF1F0H1Vertex_H
#define HERWIG_UEDF1F0H1Vertex_H
//
// This is the declaration of the UEDF1F0H1Vertex class.
//

#include "ThePEG/Helicity/Vertex/Scalar/FFSVertex.h"
#include "UEDBase.h"

namespace Herwig {
using namespace ThePEG;

/**
 * The implementation of the \f$ f_{(1)}f_{(0)}H_{(1)}^\pm\f$ vertex. It inherits 
 * from FFSVertex and implements the setCoupling member function. 
 *
 * @see \ref UEDF1F0H1VertexInterfaces "The interfaces"
 * defined for UEDF1F0H1Vertex.
 */
class UEDF1F0H1Vertex: public FFSVertex {

public:

  /**
   * The default constructor.
   */
  UEDF1F0H1Vertex();

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
  static ClassDescription<UEDF1F0H1Vertex> initUEDF1F0H1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDF1F0H1Vertex & operator=(const UEDF1F0H1Vertex &);

private:

  /**
   * The compactification radius.
   */
  InvEnergy theRadius;

  /**
   * The mass of the \f$W\f$-boson.
   */
  Energy theMw;

  /**
   * The mass of the \f$Z\f$-boson.
   */
  Energy theMz;

  /**
   * The value of \f$\sin\Theta_W\f$.
   */
  double theSinThetaW;

  /**
   * The value of \f$\cos\Theta_W\f$.
   */
  double theCosThetaW;

  /**
   * The scale at which the coupling was last evaluated. 
   */
  Energy2 theq2Last;
  
  /**
   * The value of the coupling when it was last evaluated.
   */
  Complex theCoupLast;

   /**
   * The value of the left coupling when it was last evaluated.
   */
  Complex theLeftLast;

  /**
   * The value of the right coupling when it was last evaluated.
   */
  Complex theRightLast;

  /**
   * The PDG code of the 1st particle at the vertex when it was last evaluated
   */
  long theAntiLast;

  /**
   * The PDG code of the 2nd particle at the vertex when it was last evaluated
   */
  long theFermLast;

  /**
   * The PDG code of the Higgs  particle at the vertex when it was last 
   * evaluated
   */
  long theHLast;
};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDF1F0H1Vertex. */
template <>
struct BaseClassTrait<Herwig::UEDF1F0H1Vertex,1> {
  /** Typedef of the first base class of UEDF1F0H1Vertex. */
  typedef ThePEG::Helicity::FFSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDF1F0H1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::UEDF1F0H1Vertex>
  : public ClassTraitsBase<Herwig::UEDF1F0H1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::UEDF1F0H1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDF1F0H1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDF1F0H1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#endif /* HERWIG_UEDF1F0H1Vertex_H */
