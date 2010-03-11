// -*- C++ -*-
#ifndef HERWIG_VJetGammaHardGenerator_H
#define HERWIG_VJetGammaHardGenerator_H
//
// This is the declaration of the VJetGammaHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The VJetGammaHardGenerator class generates a hard photon from
 * a vector boson plus jet event as parto f othe correction to
 * vector boson+gamma events.
 *
 * @see \ref VJetGammaHardGeneratorInterfaces "The interfaces"
 * defined for VJetGammaHardGenerator.
 */
class VJetGammaHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  VJetGammaHardGenerator();

  /**
   *  Members which must be overridden in the inheriting classes
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
  static ClassDescription<VJetGammaHardGenerator> initVJetGammaHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VJetGammaHardGenerator & operator=(const VJetGammaHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaEM_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pTmin_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VJetGammaHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VJetGammaHardGenerator,1> {
  /** Typedef of the first base class of VJetGammaHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VJetGammaHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VJetGammaHardGenerator>
  : public ClassTraitsBase<Herwig::VJetGammaHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VJetGammaHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VJetGammaHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class VJetGammaHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VJetGammaHardGenerator_H */
