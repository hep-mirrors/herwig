// -*- C++ -*-
#ifndef HERWIG_QTildeSudakovIntegrator_H
#define HERWIG_QTildeSudakovIntegrator_H
//
// This is the declaration of the QTildeSudakovIntegrator class.
//

#include "ThePEG/Pointer/ReferenceCounted.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Utilities/GaussianIntegrator.h"


namespace Herwig {

using namespace ThePEG;

class QTildeSudakovIntegrator; 

/**
 *  Struct for the inner integrand of the Sudakov
 */
struct InnerSudakovIntegrand {
  
  /**
   *  Constructor
   */
  inline InnerSudakovIntegrand(Ptr<QTildeSudakovIntegrator>::pointer iin)
    : integrator(iin) {}

  /**
   *  Default Constructor
   */
  inline InnerSudakovIntegrand() {}
  
  /**
   *  Pointer to the Integrator
   */
  Ptr<QTildeSudakovIntegrator>::pointer integrator;
  
  /**
   *  The integrand
   */
  inline double operator() (double z) const;
  /** Return type for the GaussianIntegrator */
  typedef double ValType;
  /** Argument type for the GaussianIntegrator */
  typedef double ArgType;
};

/**
 * The QTildeSudakovIntegrator class provides a method for integrating the
 * Sudakov
 */
class QTildeSudakovIntegrator: public Pointer::ReferenceCounted {

public:

  /**
   * The default constructor.
   */
  inline QTildeSudakovIntegrator() {}

  /**
   *  Constructor from a branching element
   */
  QTildeSudakovIntegrator(const BranchingElement &, Energy MergeScale, 
			  unsigned int jetMeasureMode);

  /**
   *  Return the value
   */
  double value(Energy qtildemax, Energy qtildemin);
  
  /**
   *  The integrand for the inner integral
   */
  double innerIntegrand(double) const;

  /**
   *  Integrand for the Sudakov form factor
   */
  double operator() (double) const;
  /** Return type for the GaussianIntegrator */
  typedef double ValType;
  /** Argument type for the GaussianIntegrator */
  typedef double ArgType;

  /**
   *  Minimum allowed scale
   */
  inline Energy minimumScale() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QTildeSudakovIntegrator & operator=(const QTildeSudakovIntegrator &);

private:

  /**
   * The coupling
   */
  ShowerAlphaPtr coupling_;
  
  /**
   * The Spltting function
   */
  SplittingFnPtr splittingFunction_;
  
  /**
   *  The particles
   */
  IdList ids_;

  /**
   *  Integrand of the inner integral
   */
  InnerSudakovIntegrand inner_;

  /**
   *  Gaussian integrator for the inner integral
   */
  GaussianIntegrator innerIntegrator_;

  /**
   *  Gaussian integrator for the outer integral
   */
  GaussianIntegrator outerIntegrator_;

  /**
   *  The masses of the particles in the current branching
   */
  vector<Energy> masses_;

  /**
   *  The mass squared of the particles in the current branching
   */
  vector<Energy2> masssquared_;

  /**
   *  The minimum \f$p_T\f$ for the branching
   */
  Energy pTmin_;

  /**
   *  Maximum scale for the Sudakov being integrated
   */
  Energy qtildeh_;

  /**
   *  Minimum scale for the Sudakov being integrated
   */
  Energy qtildemin_;
  
  /**
   *  The value of \f$\tilde{q}\f$ for the current integral
   */
  mutable Energy qtilde_;
  
  /**
   *  The CKKW merge scale (durham kt)
   */
  mutable Energy mergeScale_;

  /**
   *  The CKKW jet definition begin used
   */
  mutable unsigned int jetMeasureMode_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QTildeSudakovIntegrator. */
template <>
struct BaseClassTrait<Herwig::QTildeSudakovIntegrator,1> {
  /** Typedef of the first base class of QTildeSudakovIntegrator. */
  typedef Pointer::ReferenceCounted NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QTildeSudakovIntegrator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QTildeSudakovIntegrator>
  : public ClassTraitsBase<Herwig::QTildeSudakovIntegrator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QTildeSudakovIntegrator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QTildeSudakovIntegrator is implemented. It may also include several, space-separated,
   * libraries if the class QTildeSudakovIntegrator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPowhegShower.so"; }
};

/** @endcond */

}

#include "QTildeSudakovIntegrator.icc"

#endif /* HERWIG_QTildeSudakovIntegrator_H */
