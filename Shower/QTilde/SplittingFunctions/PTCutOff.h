// -*- C++ -*-
#ifndef Herwig_PTCutOff_H
#define Herwig_PTCutOff_H
//
// This is the declaration of the PTCutOff class.
//

#include "SudakovCutOff.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the PTCutOff class.
 *
 * @see \ref PTCutOffInterfaces "The interfaces"
 * defined for PTCutOff.
 */
class PTCutOff: public SudakovCutOff {

public:
  
  /**
   * The default constructor.
   */
  PTCutOff() {}

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

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
  PTCutOff & operator=(const PTCutOff &);

};

}

#endif /* Herwig_PTCutOff_H */
