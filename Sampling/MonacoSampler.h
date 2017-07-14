// -*- C++ -*-
//
// MonacoSampler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MonacoSampler_H
#define Herwig_MonacoSampler_H
//
// This is the declaration of the MonacoSampler class.
//

#include "Herwig/Sampling/BinSampler.h"

#include "Herwig/Utilities/XML/Element.h"

// work around a Boost 1.64 bug where ublas headers would fail otherwise
#include <boost/version.hpp>
#if (BOOST_VERSION / 100 >= 1064)
#include <boost/serialization/array_wrapper.hpp>
#endif

#include <boost/numeric/ublas/matrix.hpp>

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Michael Rauch
 *
 * \brief MonacoSampler samples XCombs bins using using Monaco
 *
 * @see \ref MonacoSamplerInterfaces "The interfaces"
 * defined for MonacoSampler.
 */
class MonacoSampler: public Herwig::BinSampler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MonacoSampler();

  /**
   * The destructor.
   */
  virtual ~MonacoSampler();
  //@}

public:

  /**
   * Clone this object.
   */
  Ptr<MonacoSampler>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MonacoSampler>::ptr>(clone());
  }

public:


  /**
   * Generate the next point and return its weight; store the point in
   * lastPoint().
   */
  virtual double generate();

  /**
   * Adapt this sampler after an iteration has been run
   */
  virtual void adapt();

  /**
   * Return true, if grid data exists for this sampler.
   */
  virtual bool existsGrid() const;

  /**
   * Save grid data
   */
  virtual void saveGrid() const;

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Finalize this sampler.
   */
  virtual void finalize(bool); 

  /**
   * Fill Monaco grid data from an XML element
   */
  virtual void fromXML(const XML::Element&);

  /**
   * Return an XML element for the data of the Monaco grid
   */
  virtual XML::Element toXML() const;


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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

   /**
    * Rate of grid modification (0 for no modification)
    */
  double theAlpha;

   /**
    * Number of grid divisions per dimension
    */
  size_t theGridDivisions;

   /**
    * Grid boundaries 
    * (first index: dimension of random numbers,
    * second index: dimension of partitions per random number) 
    */
  boost::numeric::ublas::matrix<double> theGrid;

   /**
    * Collected value per grid bin
    */
  boost::numeric::ublas::matrix<double> theGridData;

private:

   /**
    * Number of points collected in iteration so far
    */
  size_t theIterationPoints;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MonacoSampler & operator=(const MonacoSampler &);

};

}

#endif /* Herwig_MonacoSampler_H */
