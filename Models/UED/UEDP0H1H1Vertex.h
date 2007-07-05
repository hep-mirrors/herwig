// -*- C++ -*-
#ifndef HERWIG_UEDP0H1H1Vertex_H
#define HERWIG_UEDP0H1H1Vertex_H
//
// This is the declaration of the UEDP0H1H1Vertex class.
//

#include "Herwig++/Helicity/Vertex/Scalar/VSSVertex.h"
#include "Herwig++/Models/UED/UEDBase.h"
#include "UEDP0H1H1Vertex.fh"

namespace Herwig {
namespace Helicity {
using namespace ThePEG;

/**
 * The implementation of the \f$A^\mu_{(0)}H_{(1)}^+H_{(1)}^-\f$ vertex. This
 * class inherits from VSSVertex and implements the setCoupling member function.
 *
 * @see \ref UEDP0H1H1VertexInterfaces "The interfaces"
 * defined for UEDP0H1H1Vertex.
 */
class UEDP0H1H1Vertex: public VSSVertex {

public:

  /**
   * The default constructor.
   */
  inline UEDP0H1H1Vertex();

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
  static ClassDescription<UEDP0H1H1Vertex> initUEDP0H1H1Vertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  UEDP0H1H1Vertex & operator=(const UEDP0H1H1Vertex &);

private:

  /**
   * A pointer to the UEDBase object .
   */
  tUEDBasePtr theUEDBase;

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
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of UEDP0H1H1Vertex. */
template <>
struct BaseClassTrait<Herwig::Helicity::UEDP0H1H1Vertex,1> {
  /** Typedef of the first base class of UEDP0H1H1Vertex. */
  typedef Herwig::Helicity::VSSVertex NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the UEDP0H1H1Vertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Helicity::UEDP0H1H1Vertex>
  : public ClassTraitsBase<Herwig::Helicity::UEDP0H1H1Vertex> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::UEDP0H1H1Vertex"; }
  /**
   * The name of a file containing the dynamic library where the class
   * UEDP0H1H1Vertex is implemented. It may also include several, space-separated,
   * libraries if the class UEDP0H1H1Vertex depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwUED.so"; }
};

/** @endcond */

}

#include "UEDP0H1H1Vertex.icc"

#endif /* HERWIG_UEDP0H1H1Vertex_H */
