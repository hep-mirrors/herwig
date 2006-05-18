// -*- C++ -*-
#ifndef HERWIG_Statistic_H
#define HERWIG_Statistic_H
//
// This is the declaration of the Statistic class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Statistic.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The Statistic class is a simple class designed to 
 * store a variable for statistical analysis
 *
 * @see \ref StatisticInterfaces "The interfaces"
 * defined for Statistic.
 */
class Statistic: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Statistic();

  /**
   * The copy constructor.
   */
  inline Statistic(const Statistic &);

  /**
   * The destructor.
   */
  virtual ~Statistic();
  //@}

public:

  /**
   *  The minimum value
   */
  inline double minimum();

  /**
   *  The maximum value
   */
  inline double maximum();

  /**
   *  Operator to add another point
   */
  inline void operator+=(double);

  /**
   *  Number of points
   */
  inline unsigned int numberOfPoints();
  
  /**
   *  Mean
   */
  inline double mean();

  /**
   *  Standard Deviation
   */
  inline double stdDev();

  /**
   *  Variance
   */
  inline double var();

  /**
   *  Total entry
   */
  inline double total();

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Statistic> initStatistic;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Statistic & operator=(const Statistic &);

private:

  /**
   *   Number of entries
   */
  unsigned int _n;

  /**
   *  Sum of the values
   */ 
  double _xsum;

  /**
   *  Sum of the squares of the values
   */
  double _x2sum;

  /**
   *  The minimum value
   */
  double _min;
  
  /**
   *  The maximum value
   */
  double _max;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Statistic. */
template <>
struct BaseClassTrait<Herwig::Statistic,1> {
  /** Typedef of the first base class of Statistic. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Statistic class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Statistic>
  : public ClassTraitsBase<Herwig::Statistic> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::Statistic"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the Statistic class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "libHwUtils.so"; }
};

/** @endcond */

}

#include "Statistic.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Statistic.tcc"
#endif

#endif /* HERWIG_Statistic_H */
