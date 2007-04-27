// -*- C++ -*-
#ifndef HERWIG_UEDF1F0Z1Vertex_H
#define HERWIG_UEDF1F0Z1Vertex_H
//
// This is the declaration of the UEDF1F0Z1Vertex class.
//

#include "Herwig++/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig++/Models/UED/UEDBase.h"
#include "UEDF1F0Z1Vertex.fh"

namespace Herwig {
 namespace Helicity {
using namespace ThePEG;

/**
 * This is the implementation of the \f$ \bar{f}^{(1)} f^0 Z^{(1)}\f$
 * vertex. It inherits from FFVVertex and implements the setCoupling virtual
 * function.
 * 
 * @see \ref UEDF1F0Z1VertexInterfaces "The interfaces"
 * defined for UEDF1F0Z1Vertex.
 */
class UEDF1F0Z1Vertex: public FFVVertex {

public:

  /**
   * The default constructor.
   */
  UEDF1F0Z1Vertex();

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<UEDF1F0Z1Vertex> initUEDF1F0Z1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDF1F0Z1Vertex & operator=(const UEDF1F0Z1Vertex &);

private:
  
  /**
   * A pointer to the UEDBase object.
   */
  tUEDBasePtr theUEDBase;

  /**
   * The value of \f$\sin\theta_W\f$.
   */
  double theSinThetaW;

  /**
   * The value of \f$\cos\theta_W\f$.
   */
  double theCosThetaW;

  /**
   * The value of \f$\sin\theta_1\f$.
   */
  double theSinThetaOne;

  /**
   * The value of \f$ \cos(\theta_W - \theta_1)\f$
   */
  double theCosWmOne;
  
  /**
   * The scale at which the coupling was last evaluated.
   */
  Energy2 theq2Last;

  /**
   * The value of the coupling when it was last evaluated
   */
  Complex theCoupLast;
  
  /**
   * The ID of the last KK particle in the vertex.
   */
  long theKKLast;

  /**
   * The value of the left coupling when it was last evaluated. 
   */
  Complex theLeftLast;
  
  /**
   * The value of the right coupling when it was last evaluated. 
   */
  Complex theRightLast;
};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDF1F0Z1Vertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::UEDF1F0Z1Vertex,1> {
  /** Typedef of the first base class of UEDF1F0Z1Vertex. */
  typedef Herwig::Helicity::FFVVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDF1F0Z1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::UEDF1F0Z1Vertex>
  : public ClassTraitsBase<Herwig::Helicity::UEDF1F0Z1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::UEDF1F0Z1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDF1F0Z1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDF1F0Z1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUEDVertex.so"; }
};

/** @endcond */

}

#include "UEDF1F0Z1Vertex.icc"

#endif /* HERWIG_UEDF1F0Z1Vertex_H */
