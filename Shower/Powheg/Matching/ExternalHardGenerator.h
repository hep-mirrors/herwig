// -*- C++ -*-
#ifndef HERWIG_ExternalHardGenerator_H
#define HERWIG_ExternalHardGenerator_H
//
// This is the declaration of the ExternalHardGenerator class.
//

#include <iostream>
#include <fstream>
#include <vector>
#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Utilities/Histogram.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "PowhegHandler.h"
#include "ThePEG/Cuts/Cuts.h"



namespace Herwig {

using namespace ThePEG;
using namespace std;
/**
 * Here is the documentation of the ExternalHardGenerator class.
 *
 * @see \ref ExternalHardGeneratorInterfaces "The interfaces"
 * defined for ExternalHardGenerator.
 */
class ExternalHardGenerator: public HardestEmissionGenerator {

public:

  /**
   * The default constructor.
   */
  inline ExternalHardGenerator()
  {}

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
  inline virtual IBPtr clone() const{
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
  //@}

private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<ExternalHardGenerator> 
  initExternalHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ExternalHardGenerator & operator=(const ExternalHardGenerator &);

private:
  
  /**
   * Pointer to the CKKW handler from which the nasonTree is obtained
   */
  PowhegHandlerPtr _CKKWh;

  /**
   * Hardest emission generator with which to fill the dead zone in min mult event
   */
  HardestEmissionGeneratorPtr _hardGen;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ExternalHardGenerator. */
template <>
struct BaseClassTrait<Herwig::ExternalHardGenerator,1> {
  /** Typedef of the first base class of ExternalHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ExternalHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ExternalHardGenerator>
  : public ClassTraitsBase<Herwig::ExternalHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ExternalHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ExternalHardGenerator is implemented. 
   * It may also include several, space-separated,
   * libraries if the class ExternalHardGenerator depends
   * on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ExternalHardGenerator_H */
