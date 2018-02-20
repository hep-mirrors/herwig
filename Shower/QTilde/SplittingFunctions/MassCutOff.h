// -*- C++ -*-
#ifndef Herwig_MassCutOff_H
#define Herwig_MassCutOff_H
//
// This is the declaration of the MassCutOff class.
//

#include "SudakovCutOff.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MassCutOff class.
 *
 * @see \ref MassCutOffInterfaces "The interfaces"
 * defined for MassCutOff.
 */
class MassCutOff: public SudakovCutOff {

public:
  
  /**
   * The default constructor.
   */
  MassCutOff() {}

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MassCutOff & operator=(const MassCutOff &);

};

}

#endif /* Herwig_MassCutOff_H */
