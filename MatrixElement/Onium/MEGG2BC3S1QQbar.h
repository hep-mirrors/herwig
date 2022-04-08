// -*- C++ -*-
#ifndef Herwig_MEGG2BC3S1QQbar_H
#define Herwig_MEGG2BC3S1QQbar_H
//
// This is the declaration of the MEGG2BC3S1QQbar class.
//

#include "GGtoBCQQbarBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEGG2BC3S1QQbar class.
 *
 * @see \ref MEGG2BC3S1QQbarInterfaces "The interfaces"
 * defined for MEGG2BC3S1QQbar.
 */
class MEGG2BC3S1QQbar: public GGtoBCQQbarBase {

public:

  /**
   * The default constructor.
   */
  MEGG2BC3S1QQbar() : GGtoBCQQbarBase(541), O1_(ZERO)
  {}

public:

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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
  MEGG2BC3S1QQbar & operator=(const MEGG2BC3S1QQbar &) = delete;

private:

  /**
   *  The colour singlet matrix element
   */
  Energy3 O1_;
  
};

}

#endif /* Herwig_MEGG2BC3S1QQbar_H */
