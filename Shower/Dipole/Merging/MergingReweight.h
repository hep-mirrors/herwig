  // -*- C++ -*-
  //
  // MergingReweight.h is a part of Herwig - A multi-purpose Monte Carlo event generator
  // Copyright (C) 2002-2019 The Herwig Collaboration
  //
  // Herwig is licenced under version 3 of the GPL, see COPYING for details.
  // Please respect the MCnet academic guidelines, see GUIDELINES for details.
  //
#ifndef Herwig_MergingReweight_H
#define Herwig_MergingReweight_H
// This is the declaration of the MergingReweight class.

#include "ThePEG/MatrixElement/ReweightBase.h"

namespace Herwig {
  
  using namespace ThePEG;

/**
 * The MergingReweight class reweights subprocesses.
 *
 * @see ReweightBase
 * 
 */
class MergingReweight: public ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MergingReweight()
    : HTPower(0),MaxPTPower(0),MaxMjjPower(0), scale(50.0*GeV),onlyColoured(true) {}
  //@}

public:

  /**
   * Return the wieght for the kinematical configuation provided by
   * the assigned XComb object (in the LastXCombInfo base class).
   */
  virtual double weight() const;

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
   * Standard Init function used to initialize the interfaces.
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
   * The weight is the minimum pt/scale to a \a power.
   */
  double HTPower,MaxPTPower,MaxMjjPower;

  /**
   * The weight is the minimum pt/\a scale to a power.
   */
  Energy scale;
  
  bool onlyColoured;



private:

  /**
   * Describe a concrete base class with persistent data.
   */
  static ClassDescription<MergingReweight> initMergingReweight;

  /**
   *  Private and non-existent assignment operator.
   */
  MergingReweight & operator=(const MergingReweight &) = delete;

};





}

#endif /* Herwig_MergingReweight_H */
