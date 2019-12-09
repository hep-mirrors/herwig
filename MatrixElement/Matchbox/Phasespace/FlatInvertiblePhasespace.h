// -*- C++ -*-
//
// FlatInvertiblePhasespace.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_FlatInvertiblePhasespace_H
#define Herwig_FlatInvertiblePhasespace_H
//
// This is the declaration of the FlatInvertiblePhasespace class.
//

#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief FlatInvertiblePhasespace implements flat, invertible phase space generation.
 *
 */
class FlatInvertiblePhasespace: public MatchboxPhasespace {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FlatInvertiblePhasespace();

  /**
   * The destructor.
   */
  virtual ~FlatInvertiblePhasespace();
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
      return 1;
    return 3*nFinal - 4;
  }

public:

  /**
   * Return true, if this phase space generator is invertible
   */
  virtual bool isInvertible() const { return true; }

  /**
   * Invert the given phase space point to the random numbers which
   * would have generated it.
   */
  virtual double invertTwoToNKinematics(const vector<Lorentz5Momentum>& momenta,
					double* r) const {
    return invertKinematics(momenta,(momenta[0]+momenta[1]).m(),r);
  }

private:

  /**
   * Solve v = (n+2) * u^(n+1) - (n+1) * u^(n+2) for u
   */
  double bisect(double v, double n, 
		double target = -16., double maxLevel = 80.) const;

  /**
   * Return rho
   */
  double rho(Energy M, Energy N, Energy m) const {
    return sqrt((sqr(M)-sqr(N+m))*(sqr(M)-sqr(N-m)))/(8.*sqr(M));
  }

  /**
   * Generate intermediate masses for a massless final state
   */
  double generateIntermediates(vector<Energy>& K,
			       const double* r) const;

  /**
   * Invert intermediate masses for a massless final state
   */
  double invertIntermediates(const vector<Energy>& K,
			     double* r) const;

  /**
   * Generate intermediate masses for a massive final state
   */
  double generateIntermediates(vector<Energy>& M,
			       const vector<Energy>& m,
			       const double* r) const;

  /**
   * Invert intermediate masses for a massive final state
   */
  double invertIntermediates(const vector<Energy>& M,
			     const vector<Energy>& m,
			     double* r) const;

  /**
   * Generate momenta in the CMS
   */
  double generateKinematics(vector<Lorentz5Momentum>& P,
			    Energy Ecm,
			    const double* r) const;

  /**
   * Invert momenta in the CMS
   */
  double invertKinematics(const vector<Lorentz5Momentum>& P,
			  Energy Ecm,
			  double* r) const;

  /** 
   * Return the appropriate phase space weight, 
   * Eq. 11 in 1308.2922
   * with the factor (2 pi)^4/(2 pi)^(3n) included 
   * and the SHat of the process divided out to have everything expressed in the units of the ThePEG conventions, i.e. 
   * without the Q^2 factor
   */ 

  long double flatWeights(int n) const;

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
  FlatInvertiblePhasespace & operator=(const FlatInvertiblePhasespace &) = delete;

};

}

#endif /* Herwig_FlatInvertiblePhasespace_H */
