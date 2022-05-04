// -*- C++ -*-
#ifndef Herwig_AnalyticOmnesFunction_H
#define Herwig_AnalyticOmnesFunction_H
//
// This is the declaration of the AnalyticOmnesFunction class.
//

#include "OmnesFunction.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the AnalyticOmnesFunction class.
 *
 * @see \ref AnalyticOmnesFunctionInterfaces "The interfaces"
 * defined for AnalyticOmnesFunction.
 */
class AnalyticOmnesFunction: public OmnesFunction {

public:
  
  /**
   * The default constructor.
   */
  AnalyticOmnesFunction() : fPi_(130.7*MeV), mRho_(0.7711*GeV),
			    rhoWidth_(0.1492*GeV), rhoConst_(0.), mPi_(ZERO),
			    localParameters_(true)
  {}

  /**
   *  Method to return the function value
   */
  virtual Complex D(Energy2 s) const;

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
  AnalyticOmnesFunction & operator=(const AnalyticOmnesFunction &) = delete;

private:

  /**
   * the pion decay constant, \f$F_\pi\f$.
   */
  Energy fPi_;
  
  /**
   * The \f$\rho\f$ mass
   */
  Energy mRho_;

  /**
   * The \f$\rho\f$ width
   */
  Energy rhoWidth_;

  /**
   * Constant for the running \f$rho\f$ width.
   */
  double rhoConst_;

  /**
   * The \f$m_\pi\f$.
   */
  Energy mPi_;

  /**
   * Use local values of the parameters.
   */
  bool localParameters_;

};

}

#endif /* Herwig_AnalyticOmnesFunction_H */
