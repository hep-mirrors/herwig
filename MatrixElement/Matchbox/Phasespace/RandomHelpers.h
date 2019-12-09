// -*- C++ -*-
//
// RandomHelpers.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_RandomHelpers_H
#define HERWIG_RandomHelpers_H

#include "ThePEG/Config/ThePEG.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Phase space generation utilities.
 */
namespace RandomHelpers {

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Small helper.
 */
inline double sign(double x) {
  return x < 0. ? -1. : 1.;
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Define the generator concept.
 */
template<class Density>
struct Generator {

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const;

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const;

  /**
   * Return the density's value
   */
  double value(double x) const;

  /**
   * Return the density's normalization
   */
  double normalization() const;

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const;

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A density expression.
 */
struct Expression {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Container base class for a general density.
 */
template<>
struct Generator<Expression> {

  /**
   * The destructor.
   */
  virtual ~Generator() {}

  /**
   * Return the lower bound of the density generated.
   */
  virtual double lower() const = 0;

  /**
   * Return the upper bound of the density generated.
   */
  virtual double upper() const = 0;

  /**
   * Return the density's value
   */
  virtual double value(double x) const = 0;

  /**
   * Return the density's normalization
   */
  virtual double normalization() const = 0;

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  virtual double operator()(double r) const = 0;

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A density container.
 */
template<class Density>
struct Container {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Container class for a general density.
 */
template<class Density>
class Generator<Container<Density> >
  : public Generator<Expression> {

  /**
   * The generator.
   */
  Generator<Density> generator;

public:

  /**
   * Construct from generator.
   */
  Generator(const Generator<Density>& gen)
    : generator(gen) {}

  /**
   * Return the lower bound of the density generated.
   */
  virtual double lower() const { return generator.lower(); }

  /**
   * Return the upper bound of the density generated.
   */
  virtual double upper() const { return generator.upper(); }

  /**
   * Return the density's value
   */
  virtual double value(double x) const { return generator.value(x); }

  /**
   * Return the density's normalization
   */
  virtual double normalization() const { return generator.normalization(); }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  virtual double operator()(double r) const { return generator(r); }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a random variable and return its weight.
 */
template<class Density>
pair<double,double> generate(const Generator<Density>& gen,
			     double r) {
  double x = gen(r);

  if ( gen.value(x) != 0. )
    return make_pair(x,gen.normalization()/gen.value(x));
  else
    return make_pair(x,0.);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Remap a density to a new interval.
 */
template<class Density>
struct Remap {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a density remapped to a new interval.
 */
template<class Density>
class Generator<Remap<Density> > {

  /**
   * The underlying generator.
   */
  Generator<Density> theGenerator;

  /**
   * The new lower bound.
   */
  double theLower;

  /**
   * The new upper bound.
   */
  double theUpper;

  /**
   * Construct from generator and new boundaries.
   */
  Generator(const Generator<Density>& gen,
	    double low, double up)
    : theGenerator(gen), theLower(low), theUpper(up) {
    if ( low >= up )
      throw std::logic_error("[Generator<Remap>] Invalid boundaries.");
  }

  /**
   * Return the generator.
   */
  const Generator<Density>& generator() const { return theGenerator; }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the density's value
   */
  double value(double y) const {
    double xm = generator().lower();
    double xp = generator().upper();
    double ym = lower();
    double yp = upper();
    double x = ((xp-xm)/(yp-ym))*y+(yp*xm-ym*xp)/(yp-ym);
    return generator().value(x);
  }

  /**
   * Return the density's normalization
   */
  double normalization() const {
    double xm = generator().lower();
    double xp = generator().upper();
    double ym = lower();
    double yp = upper();
    return ((yp-ym)/(xp-xm))*generator().normalization();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const {
    double xm = generator().lower();
    double xp = generator().upper();
    double ym = lower();
    double yp = upper();
    double x = ((yp-ym)/(xp-xm))*generator()(r)+(xp*ym-xm*yp)/(xp-xm);
    return x;
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Indicate remapping of a density.
 */
struct on {

  /**
   * The new lower boundary.
   */
  double lower;

  /**
   * The new upper boundary.
   */
  double upper;

  /**
   * Construct from boundaries.
   */
  on(double a, double b)
    :lower(a), upper(b) {}
};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct a remapped density generator.
 */
template<class Density>
Generator<Remap<Density> > operator*(const Generator<Density>& gen,
				     const on& interval) {
  return Generator<Remap<Density> >(gen,interval.lower,interval.upper);
}


/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Rescale a density.
 */
template<class Density>
struct Rescale {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a rescaled density.
 */
template<class Density>
class Generator<Rescale<Density> > {

  /**
   * The underlying generator.
   */
  Generator<Density> theGenerator;

  /**
   * The rescaling factor.
   */
  double theScale;

public:

  /**
   * Construct from generator and scale.
   */
  Generator(const Generator<Density>& gen,
	    double sc)
    : theGenerator(gen), theScale(sc) {
  }

  /**
   * Return the generator.
   */
  const Generator<Density>& generator() const { return theGenerator; }

  /**
   * Return the scale
   */
  double scale() const { return theScale; }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return generator().lower(); }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return generator().upper(); }

  /**
   * Return the density's value
   */
  double value(double x) const {
    return scale()*generator().value(x);
  }

  /**
   * Return the density's normalization
   */
  double normalization() const {
    return scale()*generator().normalization();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const {
    return generator()(r);
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct a rescaled density.
 */
template<class Density>
Generator<Rescale<Density> > operator*(double a, const Generator<Density>& gen) {
  return Generator<Rescale<Density> >(gen,a);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Add two densities.
 */
template<class Density1,
	 class Density2>
struct Sum {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate the sum of two densities.
 */
template<class Density1,
	 class Density2>
class Generator<Sum<Density1,Density2> > {

  /**
   * The first generator.
   */
  Generator<Density1> theFirstGenerator;

  /**
   * The second generator.
   */
  Generator<Density2> theSecondGenerator;

  /**
   * The lower boundary.
   */
  double theLower;

  /**
   * The upper boundary.
   */
  double theUpper;

  /**
   * The fraction of the unit interval considered for the first
   * generator.
   */
  double theFraction;

public:

  /**
   * Construct from generators.
   */
  Generator(const Generator<Density1>& firstGen,
	    const Generator<Density2>& secondGen)
    : theFirstGenerator(firstGen), theSecondGenerator(secondGen),
      theLower(min(firstGen.lower(),secondGen.lower())), 
      theUpper(max(firstGen.upper(),secondGen.upper())), 
      theFraction(1.) {
    theFraction = 
      firstGenerator().normalization() / normalization();
  }

  /**
   * Return the first generator.
   */
  const Generator<Density1>& firstGenerator() const { return theFirstGenerator; }

  /**
   * Return the second generator.
   */
  const Generator<Density2>& secondGenerator() const { return theSecondGenerator; }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the fraction of the unit interval considered for the first
   * generator.
   */
  double fraction() const { return theFraction; }

  /**
   * Return the density's value
   */
  double value(double x) const {
    double res = 0.;
    if ( firstGenerator().lower() <= x &&
	 x <= firstGenerator().upper() )
      res += firstGenerator().value(x);
    if ( secondGenerator().lower() <= x &&
	 x <= secondGenerator().upper() )
      res += secondGenerator().value(x);
    return res;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const {
    return 
      firstGenerator().normalization() + secondGenerator().normalization();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const {
    return
      r < fraction() ? 
      firstGenerator()(r/fraction()) : 
      secondGenerator()((r-fraction())/(1.-fraction()));
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct the sum of two densities.
 */
template<class Density1,
	 class Density2>
Generator<Sum<Density1,Density2> > operator+(const Generator<Density1>& first,
					     const Generator<Density2>& second) {
  return Generator<Sum<Density1,Density2> >(first,second);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Indicate that the argument density should be matched to the
 * previous one in a piecewise definition.
 */
template<class Density>
struct matcher {
 
  /**
   * The generator to be matched.
   */
  Generator<Density> generator;

  /**
   * Construct from generator.
   */
  matcher(const Generator<Density>& gen)
    : generator(gen) {}

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Indicate that the argument density should be matched to the
 * previous one in a piecewise definition.
 */
template<class Density>
matcher<Density> match(const Generator<Density>& gen) {
  return matcher<Density>(gen);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct the sum of two densities, matching the first
 * summand at its upper bound to the second at its lower bound.
 */
template<class Density1,
	 class Density2>
Generator<Sum<Density1,Rescale<Density2> > > operator+(const Generator<Density1>& first,
						       const matcher<Density2>& second) {
  double matching = 
    first.value(first.upper())/
    second.generator.value(second.generator.lower());
  return Generator<Sum<Density1,Rescale<Density2> > >(first,matching*second.generator);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A piecewise defined density.
 */
template<class Density1,
	 class Density2>
struct Piecewise {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Placeholder when constructing piecewise defined densities.
 */
struct ToBeDefined {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a piecewise defined density.
 */
template<class Density1,
	 class Density2>
class Generator<Piecewise<Density1,Density2> > {

  /**
   * The first generator.
   */
  Generator<Density1> theFirstGenerator;

  /**
   * The second generator.
   */
  Generator<Density2> theSecondGenerator;

  /**
   * The lower boundary.
   */
  double theLower;

  /**
   * The transition value.
   */
  double theIntermediate;

  /**
   * The upper boundary.
   */
  double theUpper;

  /**
   * The fraction of the unit interval considered for the first
   * generator.
   */
  double theFraction;

public:

  /**
   * Construct from generators.
   */
  Generator(const Generator<Density1>& firstGen,
	    const Generator<Density2>& secondGen)
    : theFirstGenerator(firstGen), theSecondGenerator(secondGen),
      theLower(firstGen.lower()), theIntermediate(firstGen.upper()), theUpper(secondGen.upper()), 
      theFraction(1.) {
    if ( firstGenerator().upper() != secondGenerator().lower() )
      throw std::logic_error("[Generator<Piecewise>] Invalid boundaries.");
    theFraction = 
      firstGenerator().normalization() / normalization();
  }

  /**
   * Return the first generator.
   */
  const Generator<Density1>& firstGenerator() const { return theFirstGenerator; }

  /**
   * Return the second generator.
   */
  const Generator<Density2>& secondGenerator() const { return theSecondGenerator; }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the transition value.
   */
  double intermediate() const { return theIntermediate; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the fraction of the unit interval considered for the first
   * generator.
   */
  double fraction() const { return theFraction; }

  /**
   * Return the density's value
   */
  double value(double x) const {
    return
      x < intermediate() ? 
      firstGenerator().value(x) : 
      secondGenerator().value(x);
  }

  /**
   * Return the density's normalization
   */
  double normalization() const {
    return 
      firstGenerator().normalization() + secondGenerator().normalization();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const {
    return
      r < fraction() ? 
      firstGenerator()(r/fraction()) : 
      secondGenerator()((r-fraction())/(1.-fraction()));
  }

  /**
   * Construct piecewise generators.
   */
  template<class Density3>
  Generator<Piecewise<Piecewise<Density1,Density2>,Density3> >
  operator,(const Generator<Density3>& thirdGenerator) {
    return 
      Generator<Piecewise<Piecewise<Density1,Density2>,Density3> >
      (*this,thirdGenerator);
  }

  /**
   * Construct piecewise generators.
   */
  template<class Density3>
  Generator<Piecewise<Piecewise<Density1,Density2>,Rescale<Density3> > >
  operator,(const matcher<Density3>& thirdGenerator) {
    double matching = 
      value(upper())/thirdGenerator.generator.value(upper());
    return
      Generator<Piecewise<Piecewise<Density1,Density2>,Rescale<Density3> > >
      (*this,matching*thirdGenerator.generator);
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a piecewise defined density.
 */
template<class Density>
struct Generator<Piecewise<Density,ToBeDefined> > {

  /**
   * The first generator.
   */
  Generator<Density> generator;

  /**
   * Construct from generator.
   */
  Generator(const Generator<Density>& gen)
    : generator(gen) {}

  /**
   * Construct piecewise generators.
   */
  template<class Density2>
  Generator<Piecewise<Density,Density2> >
  operator,(const Generator<Density2>& secondGen) {
    return 
      Generator<Piecewise<Density,Density2> >
      (generator,secondGen);
  }

  /**
   * Construct piecewise generators.
   */
  template<class Density2>
  Generator<Piecewise<Density,Rescale<Density2> > >
  operator,(const matcher<Density2>& secondGen) {
    double matching = 
      generator.value(generator.upper())/secondGen.generator.value(generator.upper());
    return
      Generator<Piecewise<Density,Rescale<Density2> > >
      (generator,matching*secondGen.generator);
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate a piecewise defined density.
 */
template<>
struct Generator<Piecewise<ToBeDefined,ToBeDefined> > {

  /**
   * Construct piecewise generators.
   */
  template<class Density>
  Generator<Piecewise<Density,ToBeDefined> >
  operator,(const Generator<Density>& gen) {
    return 
      Generator<Piecewise<Density,ToBeDefined> >(gen);
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct a piecewise defined density.
 */
inline Generator<Piecewise<ToBeDefined,ToBeDefined> >
piecewise() {
  return Generator<Piecewise<ToBeDefined,ToBeDefined> >();
}



/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A constant density.
 */
struct Flat {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate x flat.
 */ 
template<>
class Generator<Flat> {

  /**
   * The lower boundary.
   */
  double theLower;

  /**
   * The upper boundary.
   */
  double theUpper;

public:

  /**
   * Construct from boundaries.
   */
  Generator(double low, double up)
    : theLower(low), theUpper(up) {}

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the density's value
   */
  double value(double x) const { 
    return x>=lower() && x<=upper() ? 1. : 0.;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const { return upper()-lower(); }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const { 
    return lower() + r*(upper()-lower());
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct a constant density.
 */
inline Generator<Flat> flat(double low, double up) {
  return Generator<Flat>(low,up);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A zero density.
 */
struct Zero {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate nothing.
 */ 
template<>
class Generator<Zero> {

  /**
   * The lower boundary.
   */
  double theLower;

  /**
   * The upper boundary.
   */
  double theUpper;

public:

  /**
   * Construct from boundaries.
   */
  Generator(double low, double up)
    : theLower(low), theUpper(up) {}

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the density's value
   */
  double value(double x) const { 
    return x>=lower() && x<=upper() ? Constants::epsilon : 0.;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const { return 0.; }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const { 
    return lower() + r*(upper()-lower());
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct a zero density.
 */
inline Generator<Zero> zero(double low, double up) {
  return Generator<Zero>(low,up);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief A density 1/|x-z|
 */
struct Inverse {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate x with density 1/|x-z|
 */
template<>
class Generator<Inverse> {

  /**
   * The position of the pole
   */
  double thePole;

  /**
   * The lower bound
   */
  double theLower;

  /**
   * The upper bound
   */
  double theUpper;

  /**
   * Scale for random numbers.
   */
  double theScale;

  /**
   * Offset for random mnumbers.
   */
  double theOffset;

public:

  /**
   * Construct from pole and boundaries.
   */
  Generator(double z,
	    double l, double u)
    : thePole(z),
      theLower(l), theUpper(u),
      theScale(z < l ? log((u-z)/(l-z)) : log((z-l)/(z-u))),
      theOffset(z < l ? log(l-z) : log(z-u)) {
    if ( z >= l && z <= u )
      throw std::logic_error("[Generator<Inverse>] Pole inside sampling interval.");
  }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the position of the pole.
   */
  double pole() const { return thePole; }

  /**
   * Return the scale for random numbers.
   */
  double scale() const { return theScale; }

  /**
   * Return the offset for random mnumbers.
   */
  double offset() const { return theOffset; }

  /**
   * Return the density's value
   */
  double value(double x) const {
    return x>=lower() && x<=upper() ? 1/abs(x-pole()) : 0.;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const { 
    return scale();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const { 
    return pole() + sign(upper()-pole())*exp(scale()*r+offset());
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct the density 1/|x-z|
 */
inline Generator<Inverse> inverse(double z,
				  double lower, double upper) {
  return Generator<Inverse>(z,lower,upper);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief The density |(x-z)|^p
 */
struct Power {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate x with density |(x-z)|^p
 */
template<>
class Generator<Power> {

  /**
   * The position of the pole
   */
  double thePole;

  /**
   * The power
   */
  double thePower;

  /**
   * The lower bound
   */
  double theLower;

  /**
   * The upper bound
   */
  double theUpper;

  /**
   * Scale for random numbers.
   */
  double theScale;

  /**
   * Offset for random mnumbers.
   */
  double theOffset;

public:

  /**
   * Construct from pole, power and boundaries.
   */
  Generator(double z, double p,
	    double l, double u)
    : thePole(z), thePower(p),
      theLower(l), theUpper(u),
      theScale(z<=l ? (pow(u-z,1.+p)-pow(l-z,1.+p))/(1.+p) : (pow(z-l,1.+p)-pow(z-u,1.+p))/(1.+p)),
      theOffset(z<=l ? pow(l-z,1.+p)/(1.+p) : pow(z-u,1.+p)/(1.+p)) {
    if ( p == -1. )
      throw std::logic_error("[Generator<Power>] Unit inverse. Consider using inverse().");
    if ( z >= l && z <= u && p < 0. )
      throw std::logic_error("[Generator<Power>] Pole inside sampling interval.");
    if ( z >= l && z <= u && p > 0. )
      throw std::logic_error("[Generator<Power>] Zero inside sampling interval.");
  }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the position of the pole.
   */
  double pole() const { return thePole; }

  /**
   * Return the power.
   */
  double power() const { return thePower; }

  /**
   * Return the scale for random numbers.
   */
  double scale() const { return theScale; }

  /**
   * Return the offset for random mnumbers.
   */
  double offset() const { return theOffset; }

  /**
   * Return the density's value
   */
  double value(double x) const {
    return x>=lower() && x<=upper() ? pow(abs(x-pole()),power()) : 0.;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const { 
    return scale();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const { 
    return pole() + sign(upper()-pole())*pow((1.+power())*(scale()*r+offset()),1./(1.+power()));
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct the density |(x-z)|^p
 */
inline Generator<Power> power(double z, double p,
			      double lower, double upper) {
  return Generator<Power>(z,p,lower,upper);
}

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief The density 1/((x-z)^2 + abs(w z))
 */
struct BreitWigner {};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Generate x with density 1/((x-z)^2 + abs(w z))
 */
template<>
class Generator<BreitWigner> {

  /**
   * The position of the pole
   */
  double thePole;

  /**
   * The width
   */
  double theWidth;

  /**
   * The lower bound
   */
  double theLower;

  /**
   * The upper bound
   */
  double theUpper;

  /**
   * Scale for random numbers.
   */
  double theScale;

  /**
   * Offset for random mnumbers.
   */
  double theOffset;

  /**
   * The square root of width times pole.
   */
  double theSqrtWidth;

public:

  /**
   * Construct from pole, width and boundaries.
   */
  Generator(double z, double w,
	    double l, double u)
    : thePole(z), theWidth(w),
      theLower(l), theUpper(u),
      theScale((atan((u-z)/sqrt(abs(w*z)))-atan((l-z)/sqrt(abs(w*z))))/sqrt(abs(w*z))),
      theOffset(atan((l-z)/sqrt(abs(w*z)))/sqrt(abs(w*z))),
      theSqrtWidth(sqrt(abs(w*z))) {
    if ( w == 0. )
      throw std::logic_error("[Generator<BreitWigner>] Zero width. Consider using power().");
  }

  /**
   * Return the lower bound of the density generated.
   */
  double lower() const { return theLower; }

  /**
   * Return the upper bound of the density generated.
   */
  double upper() const { return theUpper; }

  /**
   * Return the position of the pole.
   */
  double pole() const { return thePole; }

  /**
   * Return the width.
   */
  double width() const { return theWidth; }

  /**
   * Return the scale for random numbers.
   */
  double scale() const { return theScale; }

  /**
   * Return the offset for random mnumbers.
   */
  double offset() const { return theOffset; }

  /**
   * The square root of width times pole.
   */
  double sqrtWidth() const { return theSqrtWidth; }

  /**
   * Return the density's value
   */
  double value(double x) const {
    return 
      x>=lower() && x<=upper() ?
      1./(sqr(x-pole())+abs(width()*pole())) : 0.;
  }

  /**
   * Return the density's normalization
   */
  double normalization() const { 
    return scale();
  }

  /**
   * Generate the return value according to the implemented density,
   * given a flat random number on the unit interval.
   */
  double operator()(double r) const { 
    double res = pole() + sqrtWidth()*tan(sqrtWidth()*(scale()*r+offset()));
    if ( res <= lower() ) return lower()*(1+std::numeric_limits<double>::epsilon());
    else if ( res >= upper() ) return upper()*(1-std::numeric_limits<double>::epsilon());
    else return res;
  }

};

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Construct the density 1/((x-z)^2 + abs(w z))
 */
inline Generator<BreitWigner> breitWigner(double z, double w,
					  double lower, double upper) {
  return Generator<BreitWigner>(z,w,lower,upper);
}

}

}

#endif // HERWIG_RandomHelpers_H
