// -*- C++ -*-
#ifndef HERWIG_HardestEmissionGenerator_H
#define HERWIG_HardestEmissionGenerator_H
//
// This is the declaration of the HardestEmissionGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/Base/ShowerTree.fh"
#include "NasonEvolver.fh"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig++/Shower/Base/Evolver.fh"
#include "HardestEmissionGenerator.fh"
#include "HardTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HardestEmissionGenerator class.
 *
 * @see \ref HardestEmissionGeneratorInterfaces "The interfaces"
 * defined for HardestEmissionGenerator.
 */
class HardestEmissionGenerator: public Interfaced {

/**
 * The NasonEvolver is a friend to set the pointer to the Evolver
 */
friend class NasonEvolver;

public:

  /**
   * The default constructor.
   */
  inline HardestEmissionGenerator();

public:

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr)=0;

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr)=0;
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

  /**
   *  Access to the Evolver
   */
  inline tEvolverPtr evolver();

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans) throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

protected:

  /**
   *  Method to set the Evolver
   */
  void setEvolver(tEvolverPtr);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<HardestEmissionGenerator> initHardestEmissionGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardestEmissionGenerator & operator=(const HardestEmissionGenerator &);

private:

  /**
   *  Pointer to the Evolver
   */
  tEvolverPtr _evolver;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HardestEmissionGenerator. */
template <>
struct BaseClassTrait<Herwig::HardestEmissionGenerator,1> {
  /** Typedef of the first base class of HardestEmissionGenerator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HardestEmissionGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HardestEmissionGenerator>
  : public ClassTraitsBase<Herwig::HardestEmissionGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HardestEmissionGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HardestEmissionGenerator is implemented. It may also include several, space-separated,
   * libraries if the class HardestEmissionGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "HardestEmissionGenerator.icc"

#endif /* HERWIG_HardestEmissionGenerator_H */
