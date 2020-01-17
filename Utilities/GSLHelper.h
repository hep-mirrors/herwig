// -*- C++ -*-
//
// GSLHelper.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GSLHelper_H
#define HERWIG_GSLHelper_H
//
// This is the declaration of the GSLHelper class.
//

namespace Herwig {

  using namespace ThePEG;

/** \ingroup Utilities This class can be used to inherit data structures
 * from, which can then be used by the GSL algorithms that need a
 * pointer to a function and don't know about Units. This class defines
 * the necessary typedefs and forces you to define the "()" operator. In
 * addition it implements the vUnit and aUnit static methods which can
 * be overwritten if the corresponding base unit is too far from the
 * actual used unit. This removes the numerical problem that arises once
 * the base unit is several orders of magnitude away from the used unit.
 */


template <typename V, typename T> 
struct GSLHelper 
{
    
public:

  /**
   *  Constructor
   */
  GSLHelper() {}

  /**
   *  Destructor
   */
  virtual ~GSLHelper() {}

  /**
   * Typedef for Agrument type
   */
  typedef T ArgType;

  /**
   *  Typedef for Value type
   */
  typedef V ValType;

  /**
   * Value type
   */
  virtual V vUnit() const {return TypeTraits<V>::baseunit();}

  /**
   * Agrument type
   */
  virtual T aUnit() const {return TypeTraits<T>::baseunit();}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GSLHelper & operator=(const GSLHelper &) = delete;

};

}

namespace {

  template <class T> struct GSLparam {

    //The function to find root for
    const T & function;
    
  };

}

#endif /* HERWIG_GSLHelper_H */
