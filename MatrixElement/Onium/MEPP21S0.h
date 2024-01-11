// -*- C++ -*-
#ifndef Herwig_MEPP21S0_H
#define Herwig_MEPP21S0_H
//
// This is the declaration of the MEPP21S0 class.
//

#include "MEPP2OniumPowheg.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP21S0 class.
 *
 * @see \ref MEPP21S0Interfaces "The interfaces"
 * defined for MEPP21S0.
 */
class MEPP21S0: public MEPP2OniumPowheg {

public:

  /**
   * The default constructor.
   */
  MEPP21S0();

protected:

  /**
   *  The leading-order matrix element
   */
  virtual Energy2 leadingOrderME2() const;

  /**
   *   Members to calculate the real emission matrix elements
   */
  //@{
  /**
   *  The matrix element for \f$gg\to H g\f$
   */
  virtual double ggME(Energy2 s, Energy2 t, Energy2 u) const;

  /**
   *  The matrix element for \f$qg\to H q\f$
   */
  virtual double qgME(Energy2 s, Energy2 t, Energy2 u) const;

  /**
   *  The matrix element for \f$qbarg\to H qbar\f$
   */
  virtual double qbargME(Energy2 s, Energy2 t, Energy2 u) const;
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
  MEPP21S0 & operator=(const MEPP21S0 &) = delete;

private :
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy3 O1_;
  
};

}

#endif /* Herwig_MEPP21S0_H */
