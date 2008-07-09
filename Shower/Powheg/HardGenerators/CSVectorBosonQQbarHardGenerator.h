// -*- C++ -*-
#ifndef HERWIG_CSVectorBosonQQbarHardGenerator_H
#define HERWIG_CSVectorBosonQQbarHardGenerator_H
//
// This is the declaration of the CSVectorBosonQQbarHardGenerator class.
//

#include "Herwig++/Shower/Powheg/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "CSVectorBosonQQbarHardGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the CSVectorBosonQQbarHardGenerator class.
 *
 * @see \ref CSVectorBosonQQbarHardGeneratorInterfaces "The interfaces"
 * defined for CSVectorBosonQQbarHardGenerator.
 */
class CSVectorBosonQQbarHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  inline CSVectorBosonQQbarHardGenerator();

  /**
   *  Members which must be overridden in the inheriting classes
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual NasonTreePtr generateHardest(ShowerTreePtr);

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

protected:

  /**
   *  Generate the values of 
   */
  void generate(Energy2 s,double theta, double phi);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<CSVectorBosonQQbarHardGenerator> 
  initCSVectorBosonQQbarHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CSVectorBosonQQbarHardGenerator & operator=(const CSVectorBosonQQbarHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  The minimum \f$p_T\f$
   */
  Energy _ptmin;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CSVectorBosonQQbarHardGenerator. */
template <>
struct BaseClassTrait<Herwig::CSVectorBosonQQbarHardGenerator,1> {
  /** Typedef of the first base class of CSVectorBosonQQbarHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CSVectorBosonQQbarHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::CSVectorBosonQQbarHardGenerator>
  : public ClassTraitsBase<Herwig::CSVectorBosonQQbarHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::CSVectorBosonQQbarHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * CSVectorBosonQQbarHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class CSVectorBosonQQbarHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwNasonShower.so"; }
};

/** @endcond */

}

#include "CSVectorBosonQQbarHardGenerator.icc"

#endif /* HERWIG_CSVectorBosonQQbarHardGenerator_H */
