// -*- C++ -*-
#ifndef HERWIG_DefaultEmissionGenerator_H
#define HERWIG_DefaultEmissionGenerator_H
//
// This is the declaration of the DefaultEmissionGenerator class.
//

#include "HardestEmissionGenerator.h"
#include "DefaultEmissionGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DefaultEmissionGenerator class.
 *
 * @see \ref DefaultEmissionGeneratorInterfaces "The interfaces"
 * defined for DefaultEmissionGenerator.
 */
class DefaultEmissionGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  inline DefaultEmissionGenerator();

  /**
   *  Implementation of virtual functions from the base class
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual void generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
  //@}

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DefaultEmissionGenerator> initDefaultEmissionGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultEmissionGenerator & operator=(const DefaultEmissionGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultEmissionGenerator. */
template <>
struct BaseClassTrait<Herwig::DefaultEmissionGenerator,1> {
  /** Typedef of the first base class of DefaultEmissionGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultEmissionGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultEmissionGenerator>
  : public ClassTraitsBase<Herwig::DefaultEmissionGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultEmissionGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultEmissionGenerator is implemented. It may also include several, space-separated,
   * libraries if the class DefaultEmissionGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "DefaultEmissionGenerator.so"; }
};

/** @endcond */

}

#include "DefaultEmissionGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultEmissionGenerator.tcc"
#endif

#endif /* HERWIG_DefaultEmissionGenerator_H */
