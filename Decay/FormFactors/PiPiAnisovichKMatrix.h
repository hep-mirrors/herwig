// -*- C++ -*-
#ifndef Herwig_PiPiAnisovichKMatrix_H
#define Herwig_PiPiAnisovichKMatrix_H
//
// This is the declaration of the PiPiAnisovichKMatrix class.
//

#include "KMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The PiPiAnisovichKMatrix class implements the K-matrix for $\pi\pi$ scattering from Eur.Phys.J.A 16 (2003) 229-258.
 *
 * @see \ref PiPiAnisovichKMatrixInterfaces "The interfaces"
 * defined for PiPiAnisovichKMatrix.
 */
class PiPiAnisovichKMatrix: public KMatrix {

public:

  /**
   * The default constructor.
   */
  PiPiAnisovichKMatrix();

  /**
   *  Compute the K-matrix for a given scale
   * @param s The scale
   * @param Whether or not to multiply by \f$\prod_i(1-s/m^2_i\f$ to regularise the poles
   */
  virtual boost::numeric::ublas::matrix<double> K(Energy2 s, bool multiplyByPoles=false) const;

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
  PiPiAnisovichKMatrix & operator=(const PiPiAnisovichKMatrix &) = delete;

private :

  /**
   *  Parameters for the \f$K\f$-matrix
   */
  //{@
  /**
   *   \f$s_0^{\rm scat}
   */
  Energy2 s0Scatt_;

  /**
   *   Constants
   */
  vector<double> f1a_;

  /**
   *  \f$s_A\f$
   */
  double sA_;

  /**
   *  \f$s_{A0}\f$
   */
  double sA0_;

  /**
   *  Pion mass
   */
  Energy mPi_;
  //@}
};

}

#endif /* Herwig_PiPiAnisovichKMatrix_H */
