// -*- C++ -*-
//
// O2AlphaS.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_O2AlphaS_H
#define HERWIG_O2AlphaS_H
//
// This is the declaration of the O2AlphaS class.
//

#include "ThePEG/StandardModel/AlphaSBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The O2AlphaS class is the implementation of the two-loop
 * \f$\alpha_S\f$ in the same way as in FORTRAN HERWIG.
 *
 *  The input value of \f$\Lambda_{\rm QCD}\f$ is in the \f$\bar{MS}\f$
 *  scheme and can either be converted to a Monte Carlo scheme, as in the 
 *  FORTRAN program, or left in the \f$\bar{MS}\f$ scheme to evaluate
 *  the running coupling 
 *
 * @see \ref O2AlphaSInterfaces "The interfaces"
 * defined for O2AlphaS.
 */
class O2AlphaS: public AlphaSBase {

public:

  /**
   * The default constructor.
   */
  O2AlphaS() : _lambdaQCD(180.*MeV), _bcoeff(6,0.0), _ccoeff(6,0.0),
	       _lambdas(7), _threshold(6), _match(6,0.0), _copt(0) {}

  /** @name Virtual functions to override those in the base class */
  //@{
  /**
   * The \f$\alpha_S\f$. Return the QCD coupling for a given \a scale
   * using the given standard model object \a sm.
   */
  virtual double value(Energy2 scale, const StandardModelBase & sm) const;

  /**
   * Return the flavour thresholds used. The returned vector contains
   * (in position <code>i</code>) the scales when the active number of
   * flavours changes from <code>i</code> to <code>i+1</code>.
   */
  virtual vector<Energy2> flavourThresholds() const;

  /**
   * Return the \f$\Lambda_{QCD}\f$ used for different numbers of
   * active flavours.
   */
  virtual vector<Energy> LambdaQCDs() const;
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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
  O2AlphaS & operator=(const O2AlphaS &) = delete;

private:

  /**
   *  The value of \f$\Lambda_{\rm QCD}\f$ 5-flavours
   *  in the \f$\bar{\rm MS}\f$ scheme
   */
  Energy _lambdaQCD;

  /**
   *  The values of the leading-order \f$\beta\f$-function coefficients
   */
  vector<double> _bcoeff;

  /**
   *  The values of the next-to-leading-order \f$\beta\f$-function coefficients
   */
  vector<double> _ccoeff;

  /**
   *  The values of \f$\Lambda_{\rm QCD}\f$ for the diffferent number of flavours
   */
  vector<Energy> _lambdas;

  /**
   *  The flavour thresholds
   */
  vector<Energy> _threshold;
  
  /**
   *  The constants for matching
   */
  vector<double> _match;
  /**
   *  Option for the coupling
   */
  unsigned int _copt;
};

}

#endif /* HERWIG_O2AlphaS_H */
