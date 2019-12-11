// -*- C++ -*-
#ifndef Herwig_VectorCurrentDecayConstructor_H
#define Herwig_VectorCurrentDecayConstructor_H
//
// This is the declaration of the VectorCurrentDecayConstructor class.
//

#include "NBodyDecayConstructorBase.h"
#include "Herwig/Decay/WeakCurrents/WeakCurrent.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The VectorCurrentDecayConstructor class constructs the decay of low mass vector bosons via the weak currents.
 *
 * @see \ref VectorCurrentDecayConstructorInterfaces "The interfaces"
 * defined for VectorCurrentDecayConstructor.
 */
class VectorCurrentDecayConstructor: public NBodyDecayConstructorBase {

public:

  /**
   * The default constructor.
   */
  VectorCurrentDecayConstructor() : massCut_(2.*GeV)
  {}
  
  /**
   * Function used to determine allowed decaymodes, to be implemented
   * in derived class.
   *@param part vector of ParticleData pointers containing particles in model
   */
  virtual void DecayList(const set<PDPtr> & part);

  /**
   * Number of outgoing lines. Required for correct ordering (do this one next-to-last)
   */
  virtual unsigned int numBodies() const { return 999; }

  /**
   *  Cut off
   */
  Energy massCut() const { return massCut_;}

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
  VectorCurrentDecayConstructor & operator=(const VectorCurrentDecayConstructor &) = delete;

private:

  /**
   * Model Pointer
   */
  Ptr<Herwig::StandardModel>::pointer model_;

  /**
   *  Cut-off on the mass difference
   */
  Energy massCut_;

  /**
   *  The current for the mode
   */
  vector<WeakCurrentPtr> current_;
};

}

#endif /* Herwig_VectorCurrentDecayConstructor_H */
