// -*- C++ -*-
#ifndef HERWIG_GammaGammaNewHardGenerator_H
#define HERWIG_GammaGammaNewHardGenerator_H
//
// This is the declaration of the GammaGammaNewHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the GammaGammaNewHardGenerator class.
 *
 * @see \ref GammaGammaNewHardGeneratorInterfaces "The interfaces"
 * defined for GammaGammaNewHardGenerator.
 */
class GammaGammaNewHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  GammaGammaNewHardGenerator();

  /**
   *  Members overiding the base class and implementing the generation
   * of the radiation
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GammaGammaNewHardGenerator> initGammaGammaNewHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGammaNewHardGenerator & operator=(const GammaGammaNewHardGenerator &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GammaGammaNewHardGenerator. */
template <>
struct BaseClassTrait<Herwig::GammaGammaNewHardGenerator,1> {
  /** Typedef of the first base class of GammaGammaNewHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GammaGammaNewHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GammaGammaNewHardGenerator>
  : public ClassTraitsBase<Herwig::GammaGammaNewHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GammaGammaNewHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GammaGammaNewHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class GammaGammaNewHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_GammaGammaNewHardGenerator_H */
