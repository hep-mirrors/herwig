// -*- C++ -*-
//
// GeneralStatictis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_GeneralStatistics_H
#define Herwig_GeneralStatistics_H
//
// This is the declaration of the GeneralStatistics class.
//

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/Utilities/XML/Element.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief General Monte Carlo statistics.
 *
 */
class GeneralStatistics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GeneralStatistics()
    : theMaxWeight(0.), theMinWeight(Constants::MaxDouble), 
      theSumWeights(0.), theSumSquaredWeights(0.), theSumAbsWeights(0.), 
      theSelectedPoints(0), theAcceptedPoints(0),
      theNanPoints(0), theAllPoints(0),
      theLastWeight(0.) {}

  /**
   * The destructor.
   */
  virtual ~GeneralStatistics();
  //@}

public:

  /**
   * Return the last calculated chi^2.
   */
  virtual double chi2() const { return 0.; }

  /**
   * Reset these statistics.
   */
  void reset() {
    *this = GeneralStatistics();
  }

public:

  /**
   * Return the last weight encountered.
   */
  double lastWeight() const { return theLastWeight; }

  /**
   * Return the maximum absolute weight
   */
  double maxWeight() const { return theMaxWeight; }

  /**
   * Return the minimum absolute weight
   */
  double minWeight() const { return theMinWeight; }

  /**
   * Set the maximum absolute weight
   */
  void maxWeight(double w) { theMaxWeight = w; }

  /**
   * Set the minimum absolute weight
   */
  void minWeight(double w) { theMinWeight = w; }

  /**
   * Return the sum of weights
   */
  double sumWeights() const { return theSumWeights; }

  /**
   * Return the sum of squared weights 
   */
  double sumSquaredWeights() const { return theSumSquaredWeights; }

  /**
   * Return the sum of absolute weights 
   */
  double sumAbsWeights() const { return theSumAbsWeights; }

  /**
   * Return the number of selected points.
   */
  unsigned long selectedPoints() const { return theSelectedPoints; }

  /**
   * Return the nnumber of accepted points.
   */
  unsigned long acceptedPoints() const { return theAcceptedPoints; }

  /**
   * Return the number of points where a nan or inf weight has been
   * encountered.
   */
  unsigned long nanPoints() const { return theNanPoints; }

  /**
   * Return the number of all points.
   */
  unsigned long allPoints() const { return theAllPoints; }

  /**
   * Return the average weight.
   */
  virtual double averageWeight() const {
    return selectedPoints() > 0 ? sumWeights()/selectedPoints() : 0.;
  }

  /**
   * Return the average absolute weight.
   */
  virtual double averageAbsWeight() const {
    return selectedPoints() > 0 ? sumAbsWeights()/selectedPoints() : 0.;
  }

  /**
   * Return the variance of weights.
   */
  double weightVariance() const {
    return 
      selectedPoints() > 1 ? 
      abs(sumSquaredWeights() - sqr(sumWeights())/selectedPoints())/(selectedPoints()-1) : 0.;
  }

  /**
   * Return the variance of absolute weights.
   */
  double absWeightVariance() const {
    return 
      selectedPoints() > 1 ? 
      abs(sumSquaredWeights() - sqr(sumAbsWeights())/selectedPoints())/(selectedPoints()-1) : 0.;
  }

  /**
   * Return the variance of the average weight.
   */
  virtual double averageWeightVariance() const {
    return selectedPoints() > 1 ? weightVariance()/selectedPoints() : 0.;
  }

  /**
   * Return the variance of the average absolute weight.
   */
  virtual double averageAbsWeightVariance() const {
    return selectedPoints() > 1 ? absWeightVariance()/selectedPoints() : 0;
  }

  /**
   * Select an event
   */
  virtual void select(double weight, bool doIntegral = true) {
    if ( isnan(weight) || isinf(weight) ) {
      theLastWeight = weight;
      theNanPoints += 1;
      theAllPoints += 1;
      return;
    }
    theLastWeight = weight;
    theMaxWeight = max(theMaxWeight,abs(weight));
    theMinWeight = min(theMinWeight,abs(weight));
    if ( !doIntegral )
      return;
    theSumWeights += weight;
    theSumSquaredWeights += sqr(weight);
    theSumAbsWeights += abs(weight);
    theSelectedPoints += 1;
    theAllPoints += 1;
  }

  /**
   * Accept an event.
   */
  virtual void accept() {
    theAcceptedPoints += 1;
  }

  /**
   * Reject an event.
   */
  virtual void reject() {
    if ( isnan(lastWeight()) || isinf(lastWeight()) ) {
      theNanPoints -= 1;
      theAllPoints -= 1;
      return;
    }
    theSumWeights -= lastWeight();
    theSumSquaredWeights -= sqr(lastWeight());
    theSumAbsWeights -= abs(lastWeight());
    theSelectedPoints -= 1;
    theAcceptedPoints -= 1;
    theAllPoints -= 1;
  }

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
  //@}

  /**
   * Fill statistics data from an XML element
   */
  void fromXML(const XML::Element&);

  /**
   * Return an XML element for the data of this statistics
   */
  XML::Element toXML() const;

private:

  /**
   * The maximum weight encountered.
   */
  double theMaxWeight;

  /**
   * The minimum weight encountered.
   */
  double theMinWeight;

  /**
   * The sum of weights.
   */
  double theSumWeights;

  /**
   * The sum of weights squared.
   */
  double theSumSquaredWeights;

  /**
   * The sum of absolute values of weights
   */
  double theSumAbsWeights;

  /**
   * The number of selected points
   */
  unsigned long theSelectedPoints;

  /**
   * The number of accepted points
   */
  unsigned long theAcceptedPoints;

  /**
   * The number of points where an nan or inf weight was encountered.
   */
  unsigned long theNanPoints;

  /**
   * The number of all points.
   */
  unsigned long theAllPoints;

  /**
   * The last weight encountered
   */
  double theLastWeight;

};

inline PersistentOStream& operator<<(PersistentOStream& os, const GeneralStatistics& s) {
  s.put(os); return os;
}

inline PersistentIStream& operator>>(PersistentIStream& is, GeneralStatistics& s) {
  s.get(is); return is;
}

}

#endif /* Herwig_GeneralStatistics_H */
