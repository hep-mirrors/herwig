// -*- C++ -*-
//
// FFMassiveInvertedTildeKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FFMassiveInvertedTildeKinematics_H
#define HERWIG_FFMassiveInvertedTildeKinematics_H
//
// This is the declaration of the FFMassiveInvertedTildeKinematics class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer, Stephen Webster
 *
 * \brief FFMassiveInvertedTildeKinematics inverts the final-final tilde
 * kinematics.
 *
 */
class FFMassiveInvertedTildeKinematics: public Herwig::InvertedTildeKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FFMassiveInvertedTildeKinematics();

  /**
   * The destructor.
   */
  virtual ~FFMassiveInvertedTildeKinematics();
  //@}

public:

  /**
   * Perform the mapping of the tilde kinematics for the
   * last selected process and store all dimensionless
   * variables in the subtractionParameters() vector.
   * Return false, if the calculation of the real
   * kinematics was impossible for the selected configuration
   * and true on success.
   */
  virtual bool doMap(const double *);

  /**
   * Return the pt associated to the last generated splitting.
   */
  virtual Energy lastPt() const;

  /**
   * Return the momentum fraction associated to the last splitting.
   */
  virtual double lastZ() const;

  /**
   * Return the upper bound on pt
   */
  virtual Energy ptMax() const;

  /**
   * Given a pt, return the boundaries on z
   * Note that allowing parton masses these bounds may be too loose
   */
  virtual pair<double,double> zBounds(Energy pt, Energy hardPt = ZERO) const;
  
  /**
   * For generated pt and z, check if this point is
   * kinematically allowed
   */
  /*virtual*/ bool ptzAllowed(pair<Energy,double> ptz, vector<double>* values ) const;

  /**
   * Generate pt and z
   */
  virtual pair<Energy,double> generatePtZ(double& jac, const double * r, vector<double>* values) const;

public:
  
  /**
   * Triangular / Kallen function
   */
  template <class T>
  inline T rootOfKallen (T a, T b, T c) const {
    return sqrt( a*a + b*b + c*c - 2.*( a*b+a*c+b*c ) ); }
  
  // TODO: remove in both
  /**
   * stolen from FFMassiveKinematics.h
   * Perform a rotation on both momenta such that the first one will
   * point along the (positive) z axis. Rotate back to the original
   * reference frame by applying rotateUz(returnedVector) to each momentum.
   */
  ThreeVector<double> rotateToZ (Lorentz5Momentum& pTarget, Lorentz5Momentum& p1) {
    ThreeVector<double> oldAxis = pTarget.vect().unit();
    double ct = oldAxis.z(); double st = sqrt( 1.-sqr(ct) ); // cos,sin(theta)
    double cp = oldAxis.x()/st; double sp = oldAxis.y()/st; // cos,sin(phi)
    pTarget.setZ( pTarget.vect().mag() ); pTarget.setX( 0.*GeV ); pTarget.setY( 0.*GeV );
    Lorentz5Momentum p1old = p1;
    p1.setX(    sp*p1old.x() -    cp*p1old.y()                );
    p1.setY( ct*cp*p1old.x() + ct*sp*p1old.y() - st*p1old.z() );
    p1.setZ( st*cp*p1old.x() + st*sp*p1old.y() + ct*p1old.z() );
    return oldAxis;
  }
  
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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FFMassiveInvertedTildeKinematics & operator=(const FFMassiveInvertedTildeKinematics &) = delete;


};

}

#endif /* HERWIG_FFMassiveInvertedTildeKinematics_H */
