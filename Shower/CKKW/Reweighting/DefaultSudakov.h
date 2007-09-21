// -*- C++ -*-
#ifndef HERWIG_DefaultSudakov_H
#define HERWIG_DefaultSudakov_H
//
// This is the declaration of the DefaultSudakov class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "DefaultSudakov.fh"

#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"

#include "DefaultReweighter.fh"

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * DefaultSudakov is the class which performs all the
 * numerics for Sudakov reweighting in a standard CKKW
 * approach.
 *
 *@author Simon Plaetzer
 *
 * @see \ref DefaultSudakovInterfaces "The interfaces"
 * defined for DefaultSudakov.
 */
class DefaultSudakov: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DefaultSudakov();

  /**
   * Construct giving an associated DefaultReweighter
   * a splitting function, and a branching the splitting
   * function is to be evaluated for.
   */
  inline explicit DefaultSudakov (DefaultReweighterPtr,
				  SplittingFnPtr,
				  const IdList&,
				  bool useMassive = true,
				  bool initial = false);

  /**
   * The destructor.
   */
  virtual ~DefaultSudakov();
  //@}

public:

  /**
   * Initialize this Sudakov.
   */
  void initialize (bool run = true);

  /**
   * Return the value
   * for evolving from q1 to q2
   */
  inline double operator () (Energy2 q1, Energy2 q2);

  /**
   * Return the interpolated value for evolving
   * from q to the IR cutoff.
   */
  double interpolate (Energy2 q);

  /**
   * Return the splitting function
   */
  inline tSplittingFnPtr splittingFunction () const;

protected:

  /**
   * Return a unique filename for interpolation points.
   */
  string dataFile () const;

  /**
   * Return the Reweighter
   */
  inline tDefaultReweighterPtr reweighter () const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

#ifdef HERWIG_DEBUG_CKKW_CHECK_SUDAKOVS

public:

  void dumpSudakovCalls();

#endif

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}



private:

  /**
   * The associated DefaultReweighter object
   */
  DefaultReweighterPtr _reweighter;

  /**
   * The splitting function to be used.
   */
  SplittingFnPtr _splittingFunction;

  /**
   * The branching to be considered
   */
  IdList _ids;

  /**
   * Wether the branching is space- or timelike
   */
  bool _initial;

  /**
   * Wether or not to use massive splitting functions.
   */
  bool _useMassive;

  /**
   * The interpolation points in GeV^2
   */
  double * _qvalues;

  /**
   * The integral values
   */
  double * _ivalues;

  /**
   * The spline interpolation used
   */
  gsl_spline *_spline;
  
  /**
   * The acceleration used
   */
  gsl_interp_accel *_acc;

  /**
   * Switch to control memory allocation
   */
  bool _data_allocated;

  /**
   * Switch to control memory allocation
   */
  bool _interpolation_allocated;

#ifdef HERWIG_DEBUG_CKKW_CHECK_SUDAKOVS
  vector<pair<Energy2,double> > _sudakov_calls;
#endif
  

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DefaultSudakov> initDefaultSudakov;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultSudakov & operator=(const DefaultSudakov &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultSudakov. */
template <>
struct BaseClassTrait<Herwig::DefaultSudakov,1> {
  /** Typedef of the first base class of DefaultSudakov. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultSudakov class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultSudakov>
  : public ClassTraitsBase<Herwig::DefaultSudakov> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultSudakov"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultSudakov is implemented. It may also include several, space-separated,
   * libraries if the class DefaultSudakov depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "DefaultSudakov.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultSudakov.tcc"
#endif

#endif /* HERWIG_DefaultSudakov_H */
