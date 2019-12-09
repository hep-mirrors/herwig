// -*- C++ -*-
//
// FlatInvertiblePhasespaceLabFrame.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_FlatInvertibleLabframePhasespace_H
#define Herwig_FlatInvertibleLabframePhasespace_H
//
// This is the declaration of the FlatInvertibleLabframePhasespace class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/FlatInvertiblePhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Michael Rauch
 *
 * \brief FlatInvertibleLabframePhasespace implements flat, invertible phase space generation in the lab frame
 *
 */
class FlatInvertibleLabframePhasespace: public FlatInvertiblePhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FlatInvertibleLabframePhasespace();

  /**
   * The destructor.
   */
  virtual ~FlatInvertibleLabframePhasespace();
  //@}

public:

  /**
   * Generate a phase space point and return its weight.
   */
  virtual double generateTwoToNKinematics(const double*,
					  vector<Lorentz5Momentum>& momenta);

  /**
   * Return the number of random numbers required to produce a given
   * multiplicity final state.
   */
  virtual int nDimPhasespace(int nFinal) const {
    if ( nFinal == 1 )
      return 3;
    return 3*nFinal - 2;
  }

public:

  /**
   * Return true, if this phasespace generator will generate incoming
   * partons itself.
   */
  virtual bool haveX1X2() const { return true; }

  /**
   * Return true, if this phase space generator expects
   * the incoming partons in their center-of-mass system
   */
  virtual bool wantCMS() const { return false; }

  /**
   * Invert the given phase space point to the random numbers which
   * would have generated it.
   */
  virtual double invertTwoToNKinematics(const vector<Lorentz5Momentum>& momenta,
					double* r) const;

private:

  /**
   * True if SHat should be generated flat in log(SHat/S),
   * false if SHat should be generated flat in SHat.
   */
  bool theLogSHat;

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
  FlatInvertibleLabframePhasespace & operator=(const FlatInvertibleLabframePhasespace &) = delete;

};

}

#endif /* Herwig_FlatInvertiblePhasespace_H */
