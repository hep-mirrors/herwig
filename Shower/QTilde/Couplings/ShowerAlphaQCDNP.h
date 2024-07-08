// -*- C++ -*-
//
// ShowerAlphaQCDNP.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_ShowerAlphaQCDNP_H
#define HERWIG_ShowerAlphaQCDNP_H
//
// This is the declaration of the ShowerAlphaQCDNP class.
//

#include "Herwig/Shower/ShowerAlpha.h"

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
 * @see \ref ShowerAlphaQCDNPInterfaces "The interfaces"
 * defined for ShowerAlphaQCDNP.
 */
class ShowerAlphaQCDNP: public ShowerAlpha {

public:

  /**
   * The default constructor.
   */
  ShowerAlphaQCDNP() :
    ShowerAlpha(), 
    _qmin(1*GeV), _asMaxNP(1.0), _nfMaxNP(15.0),
    _thresholds(4), _lambda(4),
    _nloop(2), _thresopt(false),
    _alphain(0.118), _tolerance(1e-10),
    _maxtry(100),
    _npPower(2.), _optInputScale(91.1876_GeV),
    _quarkBranching(false) {}

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
   * Virtual method that is supposed to return the 
   * running alpha value evaluated at the input scale.
   * @param scale The scale
   * @return The coupling
   */
  virtual double showerValue(const Energy2 scale) const {
    if ( !_quarkBranching ) {
      return value(scale);
    }
    return nfNP(scale)*value(scale);
  }

  /**
   * It returns the running coupling value evaluated at the input scale
   * multiplied by the scale factor scaleFactor().
   */
  virtual double overestimateValue() const {
    return _asMaxNP;
  }
  
  /**
   * Virtual method, which 
   * should be overridden in a derived class to provide an 
   * overestimate approximation of the alpha value. 
   */
  virtual double showerOverestimateValue() const {
    if ( !_quarkBranching ) {
      return overestimateValue();
    }
    return overestimateNfNP()*overestimateValue();
  }

  /**
   *  Return the ratio of the coupling at the scale to the overestimated value
   */
  virtual double ratio(const Energy2 scale,double factor =1.) const {
    Energy2 q2 = sqr(factor)*scale;
    return value(q2)/overestimateValue();
  }

  /**
   *  Virtual method which returns the ratio of the running alpha
   * value at the input scale to the overestimated value.
   * @param scale The scale
   * @return The ratio
   */
  virtual double showerRatio(const Energy2 scale,double factor=1.) const {

      if ( !_quarkBranching ) {
	return ratio(scale,factor);
      }
      return ratioNfNP(scale,factor)*ratio(scale,factor);
  }

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

  /**
   * Match thresholds and write alpha_s
   * specified file; arguments are
   * Q_low/GeV Q_high/GeV n_steps filename
  */
  string checkNf(string args);

   string checkscaleNP(string args);

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

  /**
   * Return the non-perturbative running flavour number, normalized to
   * the number of active flavours.
   */
  double nfNP(const Energy2 scale) const;

  /**
   * Return the modification of the scale argument
   */
  Energy scaleNP(const Energy scale) const;

  /**
   * Return the the derivative of scaleNP
   */
  double scaleNPDerivative(const Energy scale) const;

  /**
   * Return the minimum of the modification of the scale argument
   */
  Energy scaleNPMin() const;


   /**
   * Return the scale q that minimizes scaleNP (i.e. scaleNP(q)=scaleNPMin()  )
   */
  Energy scaleThatmMinimizesScaleNP() const;

  /**
   * Return the maximum of the non-perturbative running flavour number
   * normalized to the number of perturbatively active flavours.
   */
  double overestimateNfNP() const {
    return _nfMaxNP;
  }

  /**
   * Return ration to the maximum of the non-perturbative running
   * flavour number normalized to the number of perturbatively active
   * flavours.
   */
  double ratioNfNP(const Energy2 scale,double factor=1.) const {
    Energy2 q2 = sqr(factor)*scale;
    return nfNP(q2)/overestimateNfNP();
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
  pair<short, Energy> getLamNf(Energy q) const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerAlphaQCDNP & operator=(const ShowerAlphaQCDNP &) = delete;

private:

  /**
   *  Transition value
   */
  Energy _qmin;

  /**
   * Maximum value of alphas
   */ 
  double _asMaxNP;

  /**
   * Maximum of non-perturbative flavour number
   */
  double _nfMaxNP;

  /**
   *  Thresholds for the different number of flavours 
   */
  vector<Energy> _thresholds;

  /**
   *  \f$\Lambda\f$ for the different number of flavours
   */
  vector<Energy> _lambda;

  /**
   * Option for the number of loops
   */
  unsigned int _nloop;

  /**
   * Option for the threshold masses
   */
  bool _thresopt;

  /**
   * Input value of \f$alpha_S(M_Z)\f$
   */
  double _alphain;

  /**
   * Tolerance for discontinuities at the thresholds
   */
  double _tolerance;

  /**
   * Maximum number of iterations for the Newton-Raphson method to converge
   */
  unsigned int _maxtry;

  
  /**
   * Power of non-perturbative modification
   */
  double _npPower;

  /**
   * Lowest cutoff scale
   */
  Energy _absoluteCutoff;

  /**
   * An optional input scale to be used for the input alphas; if zero MZ will
   * be used out of the particle data object.
   */
  Energy _optInputScale;

  /**
   * Flag if factor for quark branching should be included
   */
  bool _quarkBranching;

  /**
   * The quark masses to be used; if not empty these masses should be
   * considered instead of the ones set in the particle data objects.
   */
  vector<Energy> _quarkMasses;
  
};

}

#endif /* HERWIG_ShowerAlphaQCDNP_H */
