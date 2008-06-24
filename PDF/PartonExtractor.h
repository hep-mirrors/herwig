// -*- C++ -*-
#ifndef HERWIG_PartonExtractor_H
#define HERWIG_PartonExtractor_H
//
// This is the declaration of the PartonExtractor class.
//

#include "ThePEG/PDF/PartonExtractor.h"
#include "PartonExtractor.fh"

namespace Herwig {

/**
 * Here is the documentation of the PartonExtractor class.
 *
 * @see \ref PartonExtractorInterfaces "The interfaces"
 * defined for PartonExtractor.
 */
class PartonExtractor: public ThePEG::PartonExtractor {

public:

  /**
   * The default constructor.
   */
  inline PartonExtractor();

public:

  /**
   * Determine the number of random numbers needed to calculate
   * \f$\hat{s}\f$ and the product of all densitiy functions.
   */
  virtual ThePEG::pair<int,int> nDims(const ThePEG::PBPair & pbins);

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(ThePEG::PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(ThePEG::PersistentIStream & is, int version);
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
  inline virtual ThePEG::IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual ThePEG::IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ThePEG::ClassDescription<PartonExtractor> initPartonExtractor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PartonExtractor & operator=(const PartonExtractor &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PartonExtractor. */
template <>
struct BaseClassTrait<Herwig::PartonExtractor,1> {
  /** Typedef of the first base class of PartonExtractor. */
  typedef ThePEG::PartonExtractor NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PartonExtractor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PartonExtractor>
  : public ClassTraitsBase<Herwig::PartonExtractor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::PartonExtractor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PartonExtractor is implemented. It may also include several, space-separated,
   * libraries if the class PartonExtractor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "PartonExtractor.so"; }
};

/** @endcond */

}

#include "PartonExtractor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartonExtractor.tcc"
#endif

#endif /* HERWIG_PartonExtractor_H */
