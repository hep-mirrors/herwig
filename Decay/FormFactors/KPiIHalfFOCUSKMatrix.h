// -*- C++ -*-
#ifndef Herwig_KPiIHalfFOCUSKMatrix_H
#define Herwig_KPiIHalfFOCUSKMatrix_H
//
// This is the declaration of the KPiIHalfFOCUSKMatrix class.
//

#include "KMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The KPiIHalfFOCUSKMatrix class implements the K-matrix fit of
 * the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the \f$I=\frac12\f$
 * component of the \f$K\pi\f$ K-matrix. 
 *
 * @see \ref KPiIHalfFOCUSKMatrixInterfaces "The interfaces"
 * defined for KPiIHalfFOCUSKMatrix.
 */
class KPiIHalfFOCUSKMatrix: public KMatrix {

public:

  /**
   * The default constructor.
   */
  KPiIHalfFOCUSKMatrix();

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
  KPiIHalfFOCUSKMatrix & operator=(const KPiIHalfFOCUSKMatrix &);

private:

  /**
   *  Constants \f$C_{ij,k}\f$ from Eqn 8
   */
  vector<double> C11_,C22_,C12_;

  /**
   *  Couplings for the resonances
   */
  vector<Energy> g_;

  /**
   *   Adler zero position
   */
  Energy2 sHalf_;

  /**
   *  Normalisation scale
   */
  Energy2 sNorm_;
};

}

#endif /* Herwig_KPiIHalfFOCUSKMatrix_H */
