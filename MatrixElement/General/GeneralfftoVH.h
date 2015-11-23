// -*- C++ -*-
#ifndef HERWIG_GeneralfftoVH_H
#define HERWIG_GeneralfftoVH_H
//
// This is the declaration of the GeneralfftoVH class.
//

#include "Herwig/MatrixElement/MEfftoVH.h"
#include "GeneralfftoVH.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GeneralfftoVH class.
 *
 * @see \ref GeneralfftoVHInterfaces "The interfaces"
 * defined for GeneralfftoVH.
 */
class GeneralfftoVH: public MEfftoVH {

public:

  /**
   *  Type of process
   */
  enum Process {Lepton,HadronWplus,HadronWminus,HadronZ};

public:

  /**
   * The default constructor.
   */
  GeneralfftoVH();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
  //@}

  /**
   *  Set up the matrix element
   */
  void setProcessInfo(Process proc, PDPtr higgs,
		      AbstractVVSVertexPtr vertex,
		      unsigned int shapeOpt);

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GeneralfftoVH> initGeneralfftoVH;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralfftoVH & operator=(const GeneralfftoVH &);

private:

  /**
   *  The vector boson
   */
  Process process_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GeneralfftoVH. */
template <>
struct BaseClassTrait<Herwig::GeneralfftoVH,1> {
  /** Typedef of the first base class of GeneralfftoVH. */
  typedef Herwig::MEfftoVH NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GeneralfftoVH class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GeneralfftoVH>
  : public ClassTraitsBase<Herwig::GeneralfftoVH> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GeneralfftoVH"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GeneralfftoVH is implemented. It may also include several, space-separated,
   * libraries if the class GeneralfftoVH depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "GeneralfftoVH.so"; }
};

/** @endcond */

}

#endif /* HERWIG_GeneralfftoVH_H */
