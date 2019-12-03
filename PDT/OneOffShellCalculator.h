// -*- C++ -*-
//
// OneOffShellCalculator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_OneOffShellCalculator_H
#define HERWIG_OneOffShellCalculator_H
//
// This is the declaration of the OneOffShellCalculator class.
//
#include "GenericMassGenerator.h"
#include "WidthCalculatorBase.h"
#include "OneOffShellCalculator.fh"
#include "Herwig/Utilities/GSLIntegrator.h"


namespace Herwig {
using namespace ThePEG;

struct OneOffShellIntegrand;

/** \ingroup PDT
 *
 *  Use another <code>WidthCalculatorBase</code> object to integrate over the
 *  mass of on of the external particles which can be off-shell for running
 *  width calculations.
 *
 * @see WidthCalculatorBase
 * @see OneOffShellIntegrand
 * 
 */
class OneOffShellCalculator: public WidthCalculatorBase {

public:

  /**
   *  The OneOffShellIntegrand is a friend to allow access to the members needed
   *  for the integration without making the members public.
   */
  friend struct OneOffShellIntegrand;

public:

  /**
   * Constructor which should be used setting all the required members.
   * @param inloc The mass which is off-shell and to be integrated over.
   * @param inwidth Pointer to the WidthGeneratorBase object which calculates
   * the partial width for a given mass of the off-shell particle.
   * @param inmass Pointer to the GenericMassGenerator for the off-shell particle.
   * @param inmin The minimum mass for the off-shell particle.
   */
  OneOffShellCalculator(int inloc,WidthCalculatorBasePtr inwidth, 
			GenericMassGeneratorPtr inmass,
			Energy inmin)
    : _themass(inloc), _minmass(inmin),_onshellwidth(inwidth),_massptr(inmass)
  {}

  /**
   * member to calculate the partial width.
   * @param scale The mass squared for the decaying particle.
   * @return The partial width.
   */
  Energy partialWidth(Energy2 scale) const;

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @param mass The new value.
   * @return The mass required.
   */
  void resetMass(int imass,Energy mass) {
    _onshellwidth->resetMass(imass,mass);
  }

  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  Energy getMass(const int imass) const {
    return _onshellwidth->getMass(imass);
  }

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  Energy otherMass(const int imass) const {
    return _onshellwidth->otherMass(imass);
  }

protected:

  /**
   * The integrand.
   * @param min The mass of the off-shell particle,
   * @return The differential rate.
   */
  Energy dGamma(Energy min) const {
    _onshellwidth->resetMass(_themass,min);
    Energy wgt = (_onshellwidth->partialWidth(_scale));
    wgt*=(_massptr->weight(min));
    return wgt;
  }

private:

  /**
   * Private and non-existent assignment operator.
   */
  OneOffShellCalculator & operator=(const OneOffShellCalculator &) = delete;

private:

  /**
   * which mass is offshell.
   */
  int _themass;

  /**
   * the minimum allowed mass.
   */
  Energy _minmass;

  /**
   * pointer to object calculating the on-shell width.
   */
  WidthCalculatorBasePtr _onshellwidth;

  /**
   * pointer to object calculating the mass of the particle.
   */
  GenericMassGeneratorPtr _massptr;

  /**
   * integrator
   */
  GSLIntegrator _integrator;

  /**
   * the mass squared of the decaying particle
   */
  mutable Energy2 _scale;

};


/** \ingroup PDT
 * Class for the integrand of a matrix element where one of the outgoing
 * particles is off-shell.This class is used by the OneOffShellCalculator class
 * to perform the integral.
 *
 * @see OneOffShellCalculator
 */
struct OneOffShellIntegrand  {

  /**
   * Constructor.
   * @param in Pointer to the OneOffShellCalculator class this is doing the 
   * integration for.
   * @param m2 The mass squared of the off-shell particle for the Jacobian 
   * transform.
   * @param mw The mass times width of the off-shell particle for the Jacobian 
   * transform.
   */
  OneOffShellIntegrand(tcOneOffShellCalculatorPtr in,Energy2 m2,Energy2 mw)
    : _integrand(in),_mass2(m2),_mwidth(mw) 
  {}

  /**
   * return the value
   */
  Energy operator ()(double x) const {
    return _integrand->dGamma(sqrt(_mass2+_mwidth*tan(x)));
  }
  /** Argument type for the GSLIntegrator */
  typedef double ArgType;
  /** Return type for the GSLIntegrator */
  typedef Energy ValType;

  /**
   * pointer to the decay integrator
   */
  tcOneOffShellCalculatorPtr _integrand;

  /**
   * The mass squared for the off-shell particle for the Jacobian transform.
   */
  Energy2 _mass2;

  /**
   * The mass times width for the off-shell particle for the Jacobian transform.
   */
  Energy2 _mwidth;
};
}

#endif /* HERWIG_OneOffShellCalculator_H */
