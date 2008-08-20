// -*- C++ -*-
//
// GSLHelper.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_GSLHelper_H
#define HERWIG_GSLHelper_H
//
// This is the declaration of the GSLHelper class.
//

#include "ThePEG/Config/PhysicalQty.h"

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
  GSLHelper() {}

  virtual ~GSLHelper() {}

  typedef T ArgType;

  typedef V ValType;

  virtual V vUnit() const {return TypeTraits<V>::baseunit;}

  virtual T aUnit() const {return TypeTraits<T>::baseunit;}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GSLHelper & operator=(const GSLHelper &);

};

}

namespace {

  template <class T> struct GSLparam {

    //The function to find root for
    const T & function;
    
  };

}

#endif /* HERWIG_GSLHelper_H */
