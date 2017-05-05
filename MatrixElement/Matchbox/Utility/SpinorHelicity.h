// -*- C++ -*-
//
// SpinorHelicity.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SpinorHelicity_H
#define HERWIG_SpinorHelicity_H

#include "ThePEG/Config/Complex.h"
#include "ThePEG/Vectors/LorentzVector.h"

#include <boost/operators.hpp>

namespace Herwig {

using namespace ThePEG;

namespace SpinorHelicity {

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Tag for |p>
   *
   */
  struct PlusSpinorTag {};

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Tag for |p]
   *
   */
  struct MinusSpinorTag {};

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Tag for <p|
   *
   */
  struct PlusConjugateSpinorTag {};

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Tag for [p|
   *
   */
  struct MinusConjugateSpinorTag {};

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Helpers for commonly encountered types.
   *
   */
  template<class Value>
  struct SpinorMultiplicationTraits {

    typedef typename BinaryOpTraits<Value,Value>::MulT ResultType;
    typedef complex<ResultType> ComplexResultType;
    typedef LorentzVector<ComplexResultType> ComplexVectorResultType;

  };

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Helpers for Weyl spinors
   *
   */
  template<class Type>
  struct WeylSpinorTraits;

  // specialize for |p>
  template<>
  struct WeylSpinorTraits<PlusSpinorTag> {

    template<class Value, class MValue>
    static pair<complex<Value>,complex<Value> > 
    components(const LorentzVector<MValue>& p) {
      if ( p.t() < ZERO ) {
	pair<complex<Value>,complex<Value> > res = 
	  components<Value,MValue>(-p);
  // do not revert to *=, breaks with XCode 5.1
	res.first = res.first * Complex(0.,1.);
	res.second = res.second * Complex(0.,1.);
	return res;
      }
      Energy pPlus = p.t() + p.x();
      if ( abs(pPlus) < 1.e-10 * GeV ) {
	return make_pair(complex<Value>(ZERO),
			 complex<Value>(sqrt(2.*p.t())));
      }
      return make_pair(complex<Value>(sqrt(pPlus)),
		       complex<Value>(p.z()/sqrt(pPlus),p.y()/sqrt(pPlus)));
    }

  };

  // specialize for |p]
  template<>
  struct WeylSpinorTraits<MinusSpinorTag> {

    template<class Value, class MValue>
    static pair<complex<Value>,complex<Value> > 
    components(const LorentzVector<MValue>& p) {
      if ( p.t() < ZERO ) {
	pair<complex<Value>,complex<Value> > res = 
	  components<Value,MValue>(-p);
  // do not revert to *=, breaks with XCode 5.1  
	res.first = res.first * Complex(0.,1.);
	res.second = res.second * Complex(0.,1.);
	return res;
      }
      Energy pPlus = p.t() + p.x();
      if ( abs(pPlus) < 1.e-10 * GeV ) {
	return make_pair(complex<Value>(sqrt(2.*p.t())),
			 complex<Value>(ZERO));
      }
      return make_pair(complex<Value>(p.z()/sqrt(pPlus),-p.y()/sqrt(pPlus)),
		       -complex<Value>(sqrt(pPlus)));
    }

  };

  // specialize for <p|
  template<>
  struct WeylSpinorTraits<PlusConjugateSpinorTag> {

    typedef PlusSpinorTag ConjugateSpinorTag;
    typedef MinusSpinorTag BarSpinorTag;

    template<class Value, class MValue>
    static pair<complex<Value>,complex<Value> > 
    components(const LorentzVector<MValue>& p) {
      pair<complex<Value>,complex<Value> > res =
	WeylSpinorTraits<PlusSpinorTag>::template components<Value>(p);
      res.first = -res.first;
      swap(res.first,res.second);
      return res;
    }

  };

  // specialize for [p|
  template<>
  struct WeylSpinorTraits<MinusConjugateSpinorTag> {

    typedef MinusSpinorTag ConjugateSpinorTag;
    typedef PlusSpinorTag BarSpinorTag;

    template<class Value, class MValue>
    static pair<complex<Value>,complex<Value> > 
    components(const LorentzVector<MValue>& p) {
      pair<complex<Value>,complex<Value> > res =
	WeylSpinorTraits<MinusSpinorTag>::template components<Value>(p);
      res.second = -res.second;
      swap(res.first,res.second);
      return res;
    }

  };

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Base class for Weyl spinors
   *
   */
  template<class Type, class Value>
  class WeylSpinor {

  public:

    typedef complex<Value> ComplexType;
    typedef pair<ComplexType,ComplexType> ComponentsType;
    typedef Type Tag;
    typedef WeylSpinorTraits<Tag> Traits;
    typedef Value ValueType;

  private:

    /**
     * The components
     */
    ComponentsType theComponents;

  public:

    /**
     * Construct from components
     */
    explicit WeylSpinor(const ComponentsType& c = ComponentsType())
      : theComponents(c) {}

    /**
     * Construct from momentum
     */
    template<class MValue>
    explicit WeylSpinor(const LorentzVector<MValue>& p)
      : theComponents(Traits::template components<Value>(p)) {}

    /**
     * Return the components.
     */
    const ComponentsType& components() const { return theComponents; }

    /**
     * Return the first component
     */
    const ComplexType& s1() const { return theComponents.first; }

    /**
     * Return the second component
     */
    const ComplexType& s2() const { return theComponents.second; }

  };

  /** Define |p> */
  typedef WeylSpinor<PlusSpinorTag,SqrtEnergy> PlusSpinor;

  /** Define |p] */
  typedef WeylSpinor<MinusSpinorTag,SqrtEnergy> MinusSpinor;

  /** Define <p| */
  typedef WeylSpinor<PlusConjugateSpinorTag,SqrtEnergy> PlusConjugateSpinor;

  /** Define [p| */
  typedef WeylSpinor<MinusConjugateSpinorTag,SqrtEnergy> MinusConjugateSpinor;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Weyl spinor product
   *
   */
  template<class Type, class Value>
  class SpinorProduct 
    : public boost::addable<SpinorProduct<Type,Value> >,
      public boost::subtractable<SpinorProduct<Type,Value> >,
      public boost::multipliable<SpinorProduct<Type,Value>, double>,
      public boost::multipliable<SpinorProduct<Type,Value>, complex<double> > {

  public:

    typedef typename SpinorMultiplicationTraits<Value>::ComplexResultType ResultType;
    typedef WeylSpinor<Type,Value> LeftSpinorType;
    typedef typename WeylSpinorTraits<Type>::ConjugateSpinorTag RightSpinorTag;
    typedef WeylSpinor<RightSpinorTag,Value> RightSpinorType;

  private:

    /**
     * The result.
     */
    ResultType theResult;

  public:

    /**
     * Construct from two spinors; note that the
     * spinor metric is included, when constructing spinors.
     * Typedefs break zero products like <p|q]
     */
    explicit SpinorProduct(const LeftSpinorType& left,
			   const RightSpinorType& right)
      : theResult(left.s1()*right.s1()+left.s2()*right.s2()) {}

    /**
     * Implicitly convert to complex value
     */
    operator ResultType() const { return theResult; }

    /**
     * Return result
     */
    ResultType eval() const { return theResult; }

  public:

    SpinorProduct& operator+= (const SpinorProduct& other) {
      theResult += other.theResult;
      return *this;
    }

    SpinorProduct& operator-= (const SpinorProduct& other) {
      theResult -= other.theResult;
      return *this;
    }

    SpinorProduct& operator*= (double x) {
      theResult *= x;
      return *this;
    }

    SpinorProduct& operator*= (complex<double> x) {
      theResult *= x;
      return *this;
    }

  };

  /** Define <pq> */
  typedef SpinorProduct<PlusConjugateSpinorTag,SqrtEnergy> PlusSpinorProduct;

  /** Define [pq] */
  typedef SpinorProduct<MinusConjugateSpinorTag,SqrtEnergy> MinusSpinorProduct;

  /**
   * \ingroup Matchbox
   * \author Simon Platzer
   *
   * \brief Weyl spinor current.
   *
   */
  template<class Type, class Value>
  class SpinorCurrent 
    : public boost::addable<SpinorCurrent<Type,Value> >,
      public boost::subtractable<SpinorCurrent<Type,Value> >,
      public boost::multipliable<SpinorCurrent<Type,Value>, double>,
      public boost::multipliable<SpinorCurrent<Type,Value>, complex<double> > {

  public:

    typedef typename SpinorMultiplicationTraits<Value>::ComplexVectorResultType ResultType;
    typedef WeylSpinor<Type,Value> LeftSpinorType;
    typedef typename WeylSpinorTraits<Type>::BarSpinorTag RightSpinorTag;
    typedef WeylSpinor<RightSpinorTag,Value> RightSpinorType;

  private:

    ResultType theResult;

    /**
     * Calculate [p|\gamma^\mu|q>
     */
    ResultType evaluate(const WeylSpinor<MinusConjugateSpinorTag,Value>& left,
			const WeylSpinor<PlusSpinorTag,Value>& right) {
      return 
	ResultType(right.s1()*left.s1()-right.s2()*left.s2(),
		   complex<double>(0.,1.)*(right.s1()*left.s2()-right.s2()*left.s1()),
		   right.s1()*left.s2()+right.s2()*left.s1(),
		   right.s1()*left.s1()+right.s2()*left.s2());
    }

    /**
     * Calculate <p|\gamma^\mu|q]
     */
    ResultType evaluate(const WeylSpinor<PlusConjugateSpinorTag,Value>& left,
			const WeylSpinor<MinusSpinorTag,Value>& right) {
      return 
	ResultType(-right.s1()*left.s1()+right.s2()*left.s2(),
		   -complex<double>(0.,1.)*(right.s1()*left.s2()-right.s2()*left.s1()),
		   -right.s1()*left.s2()-right.s2()*left.s1(),
		   right.s1()*left.s1()+right.s2()*left.s2());
    }

  public:

    /**
     * Construct from two spinors.
     * Typedefs break zero products like <p|\gamma^\mu|q>
     */
    explicit SpinorCurrent(const LeftSpinorType& left,
			   const RightSpinorType& right)
      : theResult(evaluate(left,right)) {}

    /**
     * Implicitly convert to complex Lorentz vector
     */
    operator ResultType() const { return theResult; }

    /**
     * Return result
     */
    ResultType eval() const { return theResult; }

  public:

    SpinorCurrent& operator+= (const SpinorCurrent& other) {
      theResult += other.theResult;
      return *this;
    }

    SpinorCurrent& operator-= (const SpinorCurrent& other) {
      theResult -= other.theResult;
      return *this;
    }

    SpinorCurrent& operator*= (double x) {
      theResult *= x;
      return *this;
    }

    SpinorCurrent& operator*= (complex<double> x) {
      theResult *= x;
      return *this;
    }

  };

  /** Define <p|\gamma^\mu|q] */
  typedef SpinorCurrent<PlusConjugateSpinorTag,SqrtEnergy> PlusSpinorCurrent;

  /** Define [p|\gamma^\mu|q> */
  typedef SpinorCurrent<MinusConjugateSpinorTag,SqrtEnergy> MinusSpinorCurrent;

  /**
   * Return |c|^2
   */
  template<class T>
  typename BinaryOpTraits<T,T>::MulT abs2(const complex<T>& x) {
    return (x*conj(x)).real();
  }

}

}

#endif // HERWIG_SpinorHelicity_H
