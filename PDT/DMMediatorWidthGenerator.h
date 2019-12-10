// -*- C++ -*-
#ifndef Herwig_DMMediatorWidthGenerator_H
#define Herwig_DMMediatorWidthGenerator_H
//
// This is the declaration of the DMMediatorWidthGenerator class.
//

#include "ThePEG/PDT/WidthGenerator.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DMMediatorWidthGenerator class.
 *
 * @see \ref DMMediatorWidthGeneratorInterfaces "The interfaces"
 * defined for DMMediatorWidthGenerator.
 */
class DMMediatorWidthGenerator: public WidthGenerator {

public:

  /**
   * The default constructor.
   */
  DMMediatorWidthGenerator();


  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this object can be used for the given particle
   * type with the given decay map.
   */
  virtual bool accept(const ParticleData &) const;

  /**
   * Given a particle type and a mass of an instance of that particle
   * type, calculate a width.
   */
  virtual Energy width(const ParticleData &, Energy m) const;

  /**
   * Return decay map for the given particle type.
   */
  virtual DecayMap rate(const ParticleData &) const;
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

  /**
   * Overloaded function from Interfaced
   */
  virtual bool preInitialize() const {
    return true;
  }

protected:
  
  /**
   * Set the branching ratio of this mode. This requires 
   * calculating a new width for the decaying particle and reweighting
   * the current branching fractions.
   * @param dm The decaymode for which to set the branching ratio
   * @param pwidth The calculated width of the mode
   */
  void setBranchingRatio(tDMPtr dm, Energy pwidth);

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DMMediatorWidthGenerator & operator=(const DMMediatorWidthGenerator &);

private:

  /**
   *  The particle
   */
  PDPtr parent_;
  
  /**
   *  Weak currents to use for the decay
   */
  vector<WeakCurrentPtr> weakCurrents_;

};

}

#endif /* Herwig_DMMediatorWidthGenerator_H */
