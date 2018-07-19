// -*- C++ -*-
#ifndef Herwig_VariableMassCutOff_H
#define Herwig_VariableMassCutOff_H
//
// This is the declaration of the VariableMassCutOff class.
//

#include "SudakovCutOff.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the VariableMassCutOff class.
 *
 * @see \ref VariableMassCutOffInterfaces "The interfaces"
 * defined for VariableMassCutOff.
 */
class VariableMassCutOff: public SudakovCutOff {

public:

  /**
   *  Calculate the virtual masses for a branchings
   */
  virtual const vector<Energy> & virtualMasses(const IdList & ids);

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
   * The virtuality cut-off on the gluon \f$Q_g=\frac{\delta-am_q}{b}\f$
   * @param scale The scale \f$\delta\f$
   * @param mq The quark mass \f$m_q\f$.
   */
  Energy kinematicCutOff(Energy scale, Energy mq) const {
  	return max((scale -a_*mq)/b_,c_);
  }


  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VariableMassCutOff & operator=(const VariableMassCutOff &) = delete;

private:

  /**
   *  Parameters for the default Herwig cut-off option, i.e. the parameters for
   *  the \f$Q_g=\max(\frac{\delta-am_q}{b},c)\f$ kinematic cut-off
   */
  //@{
  /**
   *  The \f$a\f$ parameter
   */
  double a_ = 0.3;

  /**
   *  The \f$b\f$ parameter
   */
  double b_ = 2.3;

  /**
   *  The \f$c\f$ parameter
   */
  Energy c_ = 0.3_GeV;

  /**
   * Kinematic cutoff used in the parton shower phase space. 
   */
  Energy kinCutoffScale_ = 2.3_GeV;
  //@} 



};

}

#endif /* Herwig_VariableMassCutOff_H */
