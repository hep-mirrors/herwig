// -*- C++ -*-
#ifndef Herwig_KPiIThreeHalfFOCUSKMatrix_H
#define Herwig_KPiIThreeHalfFOCUSKMatrix_H
//
// This is the declaration of the KPiIThreeHalfFOCUSKMatrix class.
//

#include "KMatrix.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The KPiIThreeHalfFOCUSKMatrix class implements the K-matrix fit of
 * the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the \f$I=\frac32\f$
 * component of the \f$K\pi\f$ K-matrix. 
 *
 * @see \ref KPiIThreeHalfFOCUSKMatrixInterfaces "The interfaces"
 * defined for KPiIThreeHalfFOCUSKMatrix.
 */
class KPiIThreeHalfFOCUSKMatrix: public KMatrix {

public:

  /**
   * The default constructor.
   */
  KPiIThreeHalfFOCUSKMatrix();

  /**
   *  Compute the K-matrix for a given scale
   * @param s The scale
   * @param Whether ot not to multiply by \f$\prod_i(1-s/m^2_i\f$ to regularise the poles
   */
  virtual boost::numeric::ublas::matrix<double> K(Energy2 s, bool multiplyByPoles=false);

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
  KPiIThreeHalfFOCUSKMatrix & operator=(const KPiIThreeHalfFOCUSKMatrix &);

private:

  /**
   *  Constants \f$D_{22,i}\f$ from Eqn 9
   */
  vector<double> D_;

  /**
   *   Adler zero position
   */
  Energy2 sThreeHalf_;

  /**
   *  Normalisation scale
   */
  Energy2 sNorm_;
};

}

#endif /* Herwig_KPiIThreeHalfFOCUSKMatrix_H */
