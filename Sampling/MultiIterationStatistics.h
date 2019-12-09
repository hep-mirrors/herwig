// -*- C++ -*-
//
// MultiIterationStatictis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MultiIterationStatistics_H
#define Herwig_MultiIterationStatistics_H
//
// This is the declaration of the MultiIterationStatistics class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "GeneralStatistics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 * \brief Monte Carlo statistics for multiple iterations
 */
class MultiIterationStatistics: public Interfaced, public Herwig::GeneralStatistics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MultiIterationStatistics();

  /**
   * The destructor.
   */
  virtual ~MultiIterationStatistics();
  //@}

public:

  /**
   * Indicate the start of a new iteration
   */
  void nextIteration() {
    iterations().push_back(GeneralStatistics(*this));
    GeneralStatistics::reset();
  }

  /**
   * Return the iterations done so far.
   */
  const vector<GeneralStatistics>& iterations() const {
    return theIterations;
  }

  /**
   * Access the iterations done so far.
   */
  vector<GeneralStatistics>& iterations() {
    return theIterations;
  }

  /**
   * Return the last calculated chi^2.
   */
  virtual double chi2() const;

  /**
   * Return the average weight.
   */
  virtual double averageWeight(bool useAll = false) const;

  /**
   * Return the average absolute weight.
   */
  virtual double averageAbsWeight(bool useAll = false) const;

  /**
   * Return the variance of the average weight.
   */
  virtual double averageWeightVariance(bool useAll = false) const;

  /**
   * Return the variance of the average absolute weight.
   */
  virtual double averageAbsWeightVariance(bool useAll = false) const;

  /**
   * Return the minimum number of events per iteration to take this iteration
   * into account when calculating the total cross section
   */
  unsigned int minIterationPoints() const { return theMinIterationPoints; }

  /**
   * Set the minimum number of events per iteration to take this iteration
   * into account when calculating the total cross section
   */
  void minIterationPoints(unsigned int n) { theMinIterationPoints = n; }

  /**
   * Return true if integrals should be combined from all iterations
   */
  bool useAllIterations() const { return theUseAllIterations; }

  /**
   * Indicate that integrals should be combined from all iterations
   */
  void doUseAllIterations(bool yes = true) { theUseAllIterations = yes; };

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void put(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void get(PersistentIStream & is);

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
   * Fill statistics data from an XML element
   */
  void fromXML(const XML::Element&);

  /**
   * Return an XML element for the data of this statistics
   */
  XML::Element toXML() const;

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

private:

  /**
   * The currently accumulated iterations.
   */
  vector<GeneralStatistics> theIterations;

  /**
   * The minimum number of events per iteration to take this iteration
   * into account when calculating the total cross section
   */
  unsigned int theMinIterationPoints;

  /**
   * True if integrals should be combined from all iterations
   */
  bool theUseAllIterations;

};

inline PersistentOStream& operator<<(PersistentOStream& os, const MultiIterationStatistics& s) {
  s.put(os); return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, MultiIterationStatistics& s) {
  s.get(is); return is;
}

}

#endif /* Herwig_MultiIterationStatistics_H */
