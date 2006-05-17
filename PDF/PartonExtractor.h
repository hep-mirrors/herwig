// -*- C++ -*-
#ifndef HERWIG_PartonExtractor_H
#define HERWIG_PartonExtractor_H
//
// This is the declaration of the PartonExtractor class.
//

#include "ThePEG/PDF/PartonExtractor.h"
#include "PartonExtractor.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The Herwig++ parton exttactor is based on that of ThePEG with minor changes to
 * handle our different approach to the remnant
 * @see PartonExtractor
 */
class PartonExtractor: public ThePEG::PartonExtractor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline PartonExtractor();

  /**
   * The destructor.
   */
  virtual ~PartonExtractor();
  //@}

public:

  /**
   * Connect the remnants with the colour lines of the extracted
   * parton.
   */
  virtual void colourConnect(tPPtr particle, tPPtr parton,
			     const tPVector & remnants) const;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<PartonExtractor> initPartonExtractor;

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
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the PartonExtractor class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwMRST.so"; }
};

/** @endcond */

}

#include "PartonExtractor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartonExtractor.tcc"
#endif

#endif /* HERWIG_PartonExtractor_H */
