// -*- C++ -*-
#ifndef HERWIG_VGammaHardGenerator_H
#define HERWIG_VGammaHardGenerator_H
//
// This is the declaration of the VGammaHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The VGammaHardGenerator class implements the generation of the
 * hard QCD radiation in \f$q\bar q \to W^\pm/Z^0\gamma$ events.
 *
 * @see \ref VGammaHardGeneratorInterfaces "The interfaces"
 * defined for VGammaHardGenerator.
 */
class VGammaHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  VGammaHardGenerator();

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
  static ClassDescription<VGammaHardGenerator> initVGammaHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VGammaHardGenerator & operator=(const VGammaHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> _beams;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> _partons;
  //@}

  /**
   *  Whether the quark is in the + or - z direction
   */
  bool _quarkplus;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VGammaHardGenerator. */
template <>
struct BaseClassTrait<Herwig::VGammaHardGenerator,1> {
  /** Typedef of the first base class of VGammaHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VGammaHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VGammaHardGenerator>
  : public ClassTraitsBase<Herwig::VGammaHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VGammaHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VGammaHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class VGammaHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_VGammaHardGenerator_H */
