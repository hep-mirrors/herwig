// -*- C++ -*-
//
// ShowerAlphaQCD.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerAlphaQCD_H
#define HERWIG_ShowerAlphaQCD_H
//
// This is the declaration of the ShowerAlphaQCD class.
//

#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *  
 *  This concrete class provides the definition of the 
 *  pure virtual function value() and overestimateValue() for the 
 *  strong coupling.
 *
 *  A  number of different options for the running of the coupling
 *  and its initial definition are supported.
 *
 * @see \ref ShowerAlphaQCDInterfaces "The interfaces"
 * defined for ShowerAlphaQCD.
 */
class ShowerAlphaQCD: public ShowerAlpha {

public:

  /**
   * The default constructor.
   */
  ShowerAlphaQCD() : ShowerAlpha(), 
		     _qmin(0.630882*GeV), _asType(1), _asMaxNP(1.0), 
		     _thresholds(4), _lambda(4),
		     _nloop(3),_lambdaopt(false),_thresopt(false),
		     _lambdain(0.208364*GeV),_alphain(0.118),_inopt(true),_tolerance(1e-10),
		     _maxtry(100),_alphamin(0.),_val0(1.), _optInputScale(ZERO) {}

public:

  /**
   *  Methods to return the coupling
   */
  //@{
  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   * @param scale The scale
   * @return The coupling
   */
  virtual double value(const Energy2 scale) const;

  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double overestimateValue() const;

  /**
   *  Return the ratio of the coupling at the scale to the overestimated value
   */
  virtual double ratio(const Energy2 scale,double factor =1.) const;

  /**
   * Initialize this coupling.
   */
  virtual void initialize() { doinit(); }

  /**
   * A command to initialize the coupling and write
   * its value at the scale given by the argument (in GeV)
   */
  string value(string);

  /**
   * Match thresholds and write alpha_s
   * specified file; arguments are
   * Q_low/GeV Q_high/GeV n_steps filename
  */
  string check(string args);

  //@}

  /**
   *  Get the value of \f$\Lambda_{\rm QCd}\f$
   *  @param nf number of flavours
   */
  Energy lambdaQCD(unsigned int nf) {
    if      (nf <= 3)        return _lambda[0];
    else if (nf==4 || nf==5) return _lambda[nf-3];
    else                     return _lambda[3];
  }

  /**
   * Return the quark masses to be used; if not empty these masses
   * should be considered instead of the ones set in the particle data
   * objects.
   */
  const vector<Energy>& quarkMasses() const { return _quarkMasses; }

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
   *  Member functions which calculate the coupling
   */
  //@{
  /**
   * The 1,2,3-loop parametrization of \f$\alpha_S\f$.
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
  double alphaS(Energy q, Energy lam, int nf) const; 

  /**
   * The derivative of \f$\alpha_S\f$ with respect to \f$\ln(Q^2/\Lambda^2)\f$
   * @param q The scale
   * @param lam \f$\Lambda_{\rm QCD}\f$
   * @param nf The number of flavours 
   */
  double derivativealphaS(Energy q, Energy lam, int nf) const; 

  /**
   * Compute the value of \f$Lambda\f$ needed to get the input value of
   * the strong coupling at the scale given for the given number of flavours
   * using the Newton-Raphson method
   * @param match The scale for the coupling
   * @param alpha The input coupling
   * @param nflav The number of flavours
   */
  Energy computeLambda(Energy match, double alpha, unsigned int nflav) const;

  /**
   * Return the value of \f$\Lambda\f$ and the number of flavours at the scale.
   * @param q The scale
   * @return The number of flavours at the scale and \f$\Lambda\f$.
   */
  pair<short, Energy> getLamNfTwoLoop(Energy q) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlphaQCD & operator=(const ShowerAlphaQCD &) = delete;

private:

  /**
   *  Minimum value of the scale
   */
  Energy _qmin;

  /**
   *  Parameter controlling the behaviour of \f$\alpha_S\f$ in the
   *  non-perturbative region.
   */ 
  int _asType;

  /**
   *  Another parameter, a possible (maximum) value of alpha in the
   *  non-perturbative region.
   */ 
  double _asMaxNP;

  /**
   *  Thresholds for the different number of flavours 
   */
  vector<Energy> _thresholds;

  /**
   *  \f$\Lambda\f$ for the different number of flavours
   */
  vector<Energy> _lambda;

  /**
   *  Option for the number of loops
   */
  unsigned int _nloop;

  /**
   *  Option for the translation between \f$\Lambda_{\bar{MS}}\f$ and
   *  \f$\Lambda_{\rm Herwig}\f$
   */
  bool _lambdaopt;

  /**
   *  Option for the threshold masses
   */
  bool _thresopt;

  /**
   *  Input value of Lambda
   */
  Energy _lambdain;

  /**
   *  Input value of \f$alpha_S(M_Z)\f$
   */
  double _alphain;

  /**
   *  Option for the calculation of Lambda from input parameters
   */
  bool _inopt;

  /**
   *  Tolerance for discontinuities at the thresholds
   */
  double _tolerance;

  /**
   *  Maximum number of iterations for the Newton-Raphson method to converge
   */
  unsigned int _maxtry;

  /**
   *  The minimum value of the coupling
   */
  double _alphamin;

  /**
   *  Value of \f$\alpha_S\f$ at the minimum scale
   */
  double _val0;

  /**
   * An optional input scale to be used for the input alphas; if zero MZ will
   * be used out of the particle data object.
   */
  Energy _optInputScale;

  /**
   * The quark masses to be used; if not empty these masses should be
   * considered instead of the ones set in the particle data objects.
   */
  vector<Energy> _quarkMasses;
};

}

#endif /* HERWIG_ShowerAlphaQCD_H */
