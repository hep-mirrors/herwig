// -*- C++ -*-
#ifndef Herwig_DMModel_H
#define Herwig_DMModel_H
//
// This is the declaration of the DMModel class.
//

#include "Herwig/Models/General/BSMModel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DMModel class is designed to implement a simple dark matter mode
 * with fermionic dark matter and a vector mediator, as described in  arXiv:1911.11147 
 *
 * @see \ref DMModelInterfaces "The interfaces"
 * defined for DMModel.
 */
class DMModel: public BSMModel {

public:

  /**
   * The default constructor.
   */
  DMModel();

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
  DMModel & operator=(const DMModel &) = delete;

};

}

#endif /* Herwig_DMModel_H */
