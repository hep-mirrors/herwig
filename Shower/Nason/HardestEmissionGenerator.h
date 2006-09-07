// -*- C++ -*-
#ifndef HERWIG_HardestEmissionGenerator_H
#define HERWIG_HardestEmissionGenerator_H
//
// This is the declaration of the HardestEmissionGenerator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/Base/ShowerTree.fh"
#include "HardestEmissionGenerator.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HardestEmissionGenerator class.
 *
 * @see \ref HardestEmissionGeneratorInterfaces "The interfaces"
 * defined for HardestEmissionGenerator.
 */
class HardestEmissionGenerator: public Interfaced {

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
  virtual void generateHardest(ShowerTreePtr)=0;

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr)=0;
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<HardestEmissionGenerator> initHardestEmissionGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardestEmissionGenerator & operator=(const HardestEmissionGenerator &);

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
  static string library() { return "HardestEmissionGenerator.so"; }
};

/** @endcond */

}

#include "HardestEmissionGenerator.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardestEmissionGenerator.tcc"
#endif

#endif /* HERWIG_HardestEmissionGenerator_H */
