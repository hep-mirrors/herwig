// -*- C++ -*-
//
// Interpolator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Interpolator_H
#define HERWIG_Interpolator_H
//
// This is the declaration of the Interpolator class.
//

#include "ThePEG/Interface/Interfaced.h"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Utilities
 *  \author Peter Richardson
 *
 *  This class implments a polynominal interpolation of a table of values, it is
 *  based on the interpolation code in FORTRAN HERWIG. 
 *
 */

template <typename ValT, typename ArgT>
class Interpolator: public Interfaced {

public:

  /**
   *  Pointer to an Interpolator
   */
  typedef typename Ptr<Interpolator<ValT,ArgT> >::pointer Ptr;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Interpolator() : _order(3), _copyx(5),_copyfun(5) {}

  /**
   * Constructor with data as vectors.
   */
  Interpolator(const vector<ValT> & f, 
	       const vector<ArgT> & x, 
	       unsigned int order) 
    : _fun(f.size(),0.0),_xval(x.size(),0.0),_order(order),
      _funit(TypeTraits<ValT>::baseunit), 
      _xunit(TypeTraits<ArgT>::baseunit),
      _copyx(order+2),_copyfun(order+2) {
    assert(_order>0);
    assert(x.size() == f.size());
    for (size_t i = 0; i < f.size(); ++i) {
      _fun [i] = f[i] / _funit;
      _xval[i] = x[i] / _xunit;
    }
  }
  //@}

  /**
   * Constructor from bare arrays
   */
  Interpolator(size_t size, 
	       const double f[], ValT funit,
	       const double x[], ArgT xunit,
	       unsigned int order) 
    : _fun(size,0.0),_xval(size,0.0),_order(order),
      _funit(funit),_xunit(xunit), _copyx(order+2),_copyfun(order+2) {
    assert(_order>0);
    for (size_t i = 0; i < size; ++i) {
      _fun [i] = f[i];
      _xval[i] = x[i];
    }
  }
  //@}

  /**
   *  Return the interpolated value
   */
  ValT operator () (ArgT) const;
  /** Return type for GaussianIntegrator */
  typedef ValT ValType;
  /** Argument type for GaussianIntegrator */
  typedef ArgT ArgType;

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
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Interpolator<ValT,ArgT> > initInterpolator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Interpolator & operator=(const Interpolator &);
  
private:
  
  /**
   * the function values.
   */
  vector<double> _fun;

  /**
   * The x values.
   */
  vector<double> _xval;

  /**
   * the order of interpolation.
   */
  unsigned int _order;

  /**
   * The Unit of the function values
   */
  ValT _funit;

  /**
   * The Unit of the argument values
   */
  ArgT _xunit;

  /**
   *  Temporary storage vector
   */
  mutable vector<double> _copyx;

  /**
   *  Temporary storage vector
   */
  mutable vector<double> _copyfun;

};

/**
 * helper function to create InterpolatorPtr easily
 * (analogous to make_pair() )
 */
template <typename ValT, typename ArgT>
inline typename Interpolator<ValT,ArgT>::Ptr
make_InterpolatorPtr(size_t size, 
		     const double f[], ValT funit,
		     const double x[], ArgT xunit,
		     unsigned int order)
{
  return new_ptr(Interpolator<ValT,ArgT>(size,
					 f,funit,
					 x,xunit,
					 order));
}

/**
 * helper function to create InterpolatorPtr easily
 * (analogous to make_pair() )
 */
template <typename ValT, typename ArgT>
inline typename Interpolator<ValT,ArgT>::Ptr
make_InterpolatorPtr(const typename std::vector<ValT> & f, 
		     const typename std::vector<ArgT> & x, 
		     unsigned int order)
{
  return new_ptr(Interpolator<ValT,ArgT>(f,x,order));
}

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Interpolator. */
template <typename ValT, typename ArgT>
struct BaseClassTrait<Herwig::Interpolator<ValT,ArgT>,1> {
  /** Typedef of the first base class of Interpolator. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Interpolator class and the shared object where it is defined. */
template <typename ValT, typename ArgT>
struct ClassTraits<Herwig::Interpolator<ValT,ArgT> >
  : public ClassTraitsBase<Herwig::Interpolator<ValT,ArgT> > {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Interpolator<"
                               + ClassTraits<ValT>::className() + ","
                               + ClassTraits<ArgT>::className() + ">"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Interpolator is implemented. It may also include several, space-separated,
   * libraries if the class Interpolator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libHwUtils.so"; }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Interpolator.tcc"
#endif

#endif /* HERWIG_Interpolator_H */
