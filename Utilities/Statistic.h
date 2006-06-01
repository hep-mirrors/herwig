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
   * The destructor.
   */
  virtual ~Statistic();
  //@}

public:

  /**
   *  The minimum value
   */
  inline double minimum() const;

  /**
   *  The maximum value
   */
  inline double maximum() const;

  /**
   *  Operator to add another point
   */
  inline void operator+=(double);

  /**
   *  Number of points
   */
  inline unsigned int numberOfPoints() const;
  
  /**
   *  Mean
   */
  inline double mean() const;

  /**
   *  Standard Deviation
   */
  inline double stdDev() const;

  /**
   *  Variance
   */
  inline double var() const;

  /**
   *  Total entry
   */
  inline double total() const;

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
