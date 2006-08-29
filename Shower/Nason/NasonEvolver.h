// -*- C++ -*-
#ifndef HERWIG_NasonEvolver_H
#define HERWIG_NasonEvolver_H
//
// This is the declaration of the NasonEvolver class.
//

#include "Herwig++/Shower/Evolver.h"
#include "NasonEvolver.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the NasonEvolver class.
 *
 * @see \ref NasonEvolverInterfaces "The interfaces"
 * defined for NasonEvolver.
 */
class NasonEvolver: public Evolver {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline NasonEvolver();

  /**
   * The copy constructor.
   */
  inline NasonEvolver(const NasonEvolver &);

  /**
   * The destructor.
   */
  virtual ~NasonEvolver();
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
  static ClassDescription<NasonEvolver> initNasonEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NasonEvolver & operator=(const NasonEvolver &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of NasonEvolver. */
template <>
struct BaseClassTrait<Herwig::NasonEvolver,1> {
  /** Typedef of the first base class of NasonEvolver. */
  typedef Herwig::Evolver NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the NasonEvolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::NasonEvolver>
  : public ClassTraitsBase<Herwig::NasonEvolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::NasonEvolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * NasonEvolver is implemented. It may also include several, space-separated,
   * libraries if the class NasonEvolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "NasonEvolver.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "NasonEvolver.tcc"
#endif

#endif /* HERWIG_NasonEvolver_H */
