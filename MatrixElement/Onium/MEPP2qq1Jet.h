// -*- C++ -*-
#ifndef Herwig_MEPP2qq1Jet_H
#define Herwig_MEPP2qq1Jet_H
//
// This is the declaration of the MEPP2qq1Jet class.
//

#include "MEPP2DiquarkJet.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2qq1Jet class implements the \f$2\to2\f$ processes for the production
 * a spin-1 diquark with identical quarks
 * 
 *
 * @see \ref MEPP2qq1JetInterfaces "The interfaces"
 * defined for MEPP2qq1Jet.
 */
class MEPP2qq1Jet: public MEPP2DiquarkJet {

public:

  /**
   * The default constructor.
   */
  MEPP2qq1Jet() : type_(0)
  {}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
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
  MEPP2qq1Jet & operator=(const MEPP2qq1Jet &) = delete;

private:
  
  /**
   *  The \f$O_1\f$ colour-singlet coefficient
   */
  Energy3 R02_;

  /**
   *  Option for the type of diquark
   */
  unsigned int type_;

};

}

#endif /* Herwig_MEPP2qq1Jet_H */
