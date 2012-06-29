// -*- C++ -*-
#ifndef HERWIG_AdditionalGaugeParticleData_H
#define HERWIG_AdditionalGaugeParticleData_H
//
// This is the declaration of the AdditionalGaugeParticleData class.
//

#include "AdditionalGaugeParticleData.fh"
#include "ThePEG/PDT/ConstituentParticleData.h"

namespace Herwig {

using namespace ThePEG;

/**
 *  The HiddenPDT class is a helper class implementing enumerations for 
 *  the quantum numbers under the new gauge group.
 */
class HiddenPDT {

public:

  enum HiddenColour {
    HiddenColourNeutral = 0,               /**< Hidden colour-singlet */
    HiddenColourFundamental      = 3,      /**< Hidden Colour-triplet */
    HiddenColourAntiFundamental = -3,      /**< Hidden Colour-anti-triplet */
    HiddenColourAdjoint         = 8        /**< Hidden Colour-octet */
  };
};

/**
 * Here is the documentation of the AdditionalGaugeParticleData class.
 *
 * @see \ref AdditionalGaugeParticleDataInterfaces "The interfaces"
 * defined for AdditionalGaugeParticleData.
 */
class AdditionalGaugeParticleData: public ConstituentParticleData {

public:

  /**
   * The default constructor.
   */
  AdditionalGaugeParticleData() : hiddenColour_(HiddenPDT::HiddenColourNeutral) {}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(long newId, string newPDGName);

  /**
   * Create a particle - anti particle pair.
   */
  static PDPair Create(long newId, string newPDGName, string newAntiPDGName);
  //@}

public:

  /**
   *  Access to the hidden colour
   */
  HiddenPDT::HiddenColour hiddenColour() const {return hiddenColour_;}

  /**
   *  Whether or not coloured under the hidden group
   */
  bool hiddenColoured() const
  {return hiddenColour_!=HiddenPDT::HiddenColourNeutral;}

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
   * Protected constructor only to be used by subclasses or by the
   * Create method.
   */
  AdditionalGaugeParticleData(long newId, string newPDGName);

  /**
   * Read setup info from a standard stream. The information that must
   * be supplied is the same as for ParticleData::readSetup with an
   * additional constituent mass (in GeV) added in the end.
   */
  virtual void readSetup(istream & is);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<AdditionalGaugeParticleData> initAdditionalGaugeParticleData;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AdditionalGaugeParticleData & operator=(const AdditionalGaugeParticleData &);

private: 

  /**
   *  Colour under the new group
   */
  HiddenPDT::HiddenColour hiddenColour_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AdditionalGaugeParticleData. */
template <>
struct BaseClassTrait<Herwig::AdditionalGaugeParticleData,1> {
  /** Typedef of the first base class of AdditionalGaugeParticleData. */
  typedef ConstituentParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AdditionalGaugeParticleData class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::AdditionalGaugeParticleData>
  : public ClassTraitsBase<Herwig::AdditionalGaugeParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::AdditionalGaugeParticleData"; }
  /**
   * The name of a file containing the dynamic library where the class
   * AdditionalGaugeParticleData is implemented. It may also include several, space-separated,
   * libraries if the class AdditionalGaugeParticleData depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwHiddenValleyModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_AdditionalGaugeParticleData_H */
