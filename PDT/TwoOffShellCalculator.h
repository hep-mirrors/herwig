// -*- C++ -*-
//
// TwoOffShellCalculator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TwoOffShellCalculator_H
#define HERWIG_TwoOffShellCalculator_H
//
// This is the declaration of the TwoOffShellCalculator class.
//
#include "WidthCalculatorBase.h"
#include "GenericMassGenerator.h"
#include "TwoOffShellCalculator.fh"
#include "OneOffShellCalculator.fh"
#include "Herwig/Utilities/GSLIntegrator.h"

namespace Herwig {
using namespace ThePEG;

struct TwoOffShellIntegrand;

/** \ingroup PDT
 *
 *  Use <code>WidthCalculatorBase</code> objects to integrate over the mass of two
 *  external particles which can be off-shell for running width calculations.
 *
 * @see WidthCalculatorBase
 * @see TwoOffShellIntegrand
 */
class TwoOffShellCalculator: public WidthCalculatorBase {

/**
 *  The TwoOffShellIntegrand class is a friend to allow access to the private
 *  members for the integration.
 */
friend struct TwoOffShellIntegrand;

public:

  /**
   * Constructor which should be used setting all the required members.
   * @param inloc The mass which is off-shell and to be integrated over.
   * @param inwidth Pointer to the WidthGeneratorBase object which calculates
   * the partial width for a given mass of the off-shell particle. This
   * should be a OneOffShellCalculator instance.
   * @param inmass Pointer to the GenericMassGenerator for the off-shell particle.
   * @param inmin1 The minimum mass for the first off-shell particle.
   * @param inmin2 The minimum mass for the second off-shell particle.
   */
  TwoOffShellCalculator(int inloc, WidthCalculatorBasePtr inwidth,
			GenericMassGeneratorPtr inmass,
			Energy inmin2,Energy inmin1)
    : _themass(inloc),_minmass(inmin2),_mother(inmin1),_oneoffwidth(inwidth),
      _massptr(inmass)   
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
    _oneoffwidth->resetMass(imass,mass);
  }
  
  /**
   * Get the mass of one of the decay products.  This must be 
   * implemented in classes inheriting from this one.
   * @param imass The mass required.
   * @return The mass required.
   */
  Energy getMass(const int imass) const {
    return _oneoffwidth->getMass(imass);
  }

  /**
   * Get the masses of all bar the one specified. Used to get the limits
   * for integration.
   * @param imass The particle not needed
   * @return The sum of the other masses.
   */
  Energy otherMass(const int imass) const {
    return _oneoffwidth->otherMass(imass);
  }

protected:

  /**
   * The integrand.
   * @param mass The mass of the second off-shell particle,
   * @return The differential rate.
   */
  Energy dGamma(Energy mass) const {
    _oneoffwidth->resetMass(_themass,mass);
    Energy wgt = (_oneoffwidth->partialWidth(_scale));
    wgt*=(_massptr->weight(mass));
    return wgt;
  }

private:

  /**
   * Private and non-existent assignment operator.
   */
  TwoOffShellCalculator & operator=(const TwoOffShellCalculator &) = delete;

private:

  /**
   * The second mass which is offshell
   */

  int _themass;
  /**
   * the minimum allowed mass
   */
  Energy _minmass;

  /**
   * sum of the masses of the other decay products
   */
  Energy _mother;

  /**
   * pointer to object calculating the width for one-off shell particle.
   */
  WidthCalculatorBasePtr _oneoffwidth;

  /**
   * pointer to object calculating the mass of the particle
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
 * Class for the integrand of a matrix element where two of the outgoing
 * particles is off-shell. This class is used by the TwoOffShellCalculator class
 * to perform the integral.
 */
struct TwoOffShellIntegrand {
  /**
   * Constructor.
   * @param in Pointer to the OneOffShellCalculator class this is doing the 
   * integration for.
   * @param m2 The mass squared of the off-shell particle for the Jacobian 
   * transform.
   * @param mw The mass times width of the off-shell particle for the Jacobian 
   * transform.
   */
  TwoOffShellIntegrand(tcTwoOffShellCalculatorPtr in,Energy2 m2,Energy2 mw)
    : _integrand(in),_mass2(m2),_mwidth(mw)
  {}

  /**
   * Retreive function value
   */
  Energy operator ()(double x) const {
    return _integrand->dGamma(sqrt(_mass2+_mwidth*tan(x)));
  }
  /** Argument type for the GSLIntegrator */
  typedef double ArgType;
  /** Return type for the GSLIntegrator */
  typedef Energy ValType;

private:

  /**
   * pointer to the decay integrator
   */
  cTwoOffShellCalculatorPtr _integrand;

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

#endif /* HERWIG_TwoOffShellCalculator_H */
