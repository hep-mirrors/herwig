// -*- C++ -*-
//
// ExSampler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ExSampler_H
#define Herwig_ExSampler_H
//
// This is the declaration of the ExSampler class.
//

#include "BinSampler.h"
#include "Herwig++/Exsample2/exsample/generator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \brief Interface to the exsample generator.
 * \author Simon Platzer
 *
 * @see \ref ExSamplerInterfaces "The interfaces"
 * defined for ExSampler.
 */
class ExSampler: public Herwig::BinSampler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ExSampler();

  /**
   * The destructor.
   */
  virtual ~ExSampler();
  //@}

public:

  /**
   * Return true, if this bin sampler produces unweighted events.
   */
  virtual bool isUnweighting() const { return true; }

  /**
   * Return true, if this sampler is in a compensating mode.
   */
  virtual bool compensating() const { return generator_.compensating(); }

  /**
   * Generate the next point; store the point in lastPoint() and its
   * weight using select(); if noMaxInfo is true, do not throw
   * NewMaximum or UpdateCrossSections exceptions.
   */
  virtual void generate(bool noMaxInfo = false);

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Finalize this sampler.
   */
  virtual void finalize(bool verbose);

public:

  /**
   * Reject an event.
   */
  virtual void reject() {
    GeneralStatistics::reject();
    generator_.reject();
  }

public:

  /**
   * Evaluate with given random numbers.
   */
  double evaluate(const vector<double>& p) const {
    double ret;
    try {
      ret = eventHandler()->dSigDR(p) / nanobarn;
    } catch (Veto&) {
      ret = 0.0;
    } catch (...) {
      throw;
    }
    return ret;
  }

  /**
   * Return the lower left and upper right
   * corners of the support of this function
   */
  pair<vector<double>,vector<double> > support() const {
    vector<double> lower(dimension(),0.);
    vector<double> upper(dimension(),1.);
    return make_pair(lower,upper);
  }

  /**
   * Indicate start of presampling
   */
  void start_presampling() { }

  /**
   * Indicate end of presampling
   */
  void stop_presampling() { }

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
   * The number of presampling points.
   */
  unsigned long presampling_points_;

  /**
   * The number of points below which the grid is frozen
   */
  unsigned long freeze_grid_;

  /**
   * The efficiency threshold
   */
  double efficiency_threshold_;

  /**
   * The gains threshold.
   */
  double gain_threshold_;

  /**
   * The generator used.
   */
  exsample::generator<ExSampler,UseRandom> generator_;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ExSampler & operator=(const ExSampler &);

};

}

#endif /* Herwig_ExSampler_H */
