// -*- C++ -*-
#ifndef Herwig_ElectroWeakReweighter_H
#define Herwig_ElectroWeakReweighter_H
//
// This is the declaration of the ElectroWeakReweighter class.
//

#include "ThePEG/MatrixElement/ReweightBase.h"
#include "EWCouplings.h"
#include "CollinearSudakov.h"
#include "SoftSudakov.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ElectroWeakReweighter class.
 *
 * @see \ref ElectroWeakReweighterInterfaces "The interfaces"
 * defined for ElectroWeakReweighter.
 */
class ElectroWeakReweighter: public ReweightBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ElectroWeakReweighter();

  /**
   * The destructor.
   */
  virtual ~ElectroWeakReweighter();
  //@}  

public:

  /**
   * Return the weight for the kinematical configuation provided by
   * the assigned XComb object (in the LastXCombInfo base class).
   */
  virtual double weight() const;

  /**
   *
   */
  static tEWCouplingsPtr coupling() {
    assert(staticEWCouplings_);
    return staticEWCouplings_;
  }

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
   *  Check the evolution for a fixed s,t,u
   */
  void testEvolution(Energy2 s,Energy2 t, Energy2 u) const;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ElectroWeakReweighter & operator=(const ElectroWeakReweighter &);

private:

  /**
   *  The Electroweak Couplings
   */
  EWCouplingsPtr EWCouplings_;

  /**
   *  The Collinear Sudakov
   */
  CollinearSudakovPtr collinearSudakov_;

  /**
   *  The Soft Sudakov
   */
  SoftSudakovPtr softSudakov_;

  /**
   *  The couplings to allow global access
   */
  static tEWCouplingsPtr staticEWCouplings_;

};

}

#endif /* Herwig_ElectroWeakReweighter_H */
