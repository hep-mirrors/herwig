// -*- C++ -*-
//
// GeneralSampler.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_GeneralSampler_H
#define Herwig_GeneralSampler_H
//
// This is the declaration of the GeneralSampler class.
//

#include "ThePEG/Handlers/SamplerBase.h"
#include "BinSampler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief A GeneralSampler class
 *
 * @see \ref GeneralSamplerInterfaces "The interfaces"
 * defined for GeneralSampler.
 */
class GeneralSampler: public SamplerBase {

public:

  /**
   * Exception thrown in case too many attempts to unweight an event.
   */
  struct MaxTryException : public ThePEG::Exception {};

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GeneralSampler();

  /**
   * The destructor.
   */
  virtual ~GeneralSampler();
  //@}

public:

  /** @name Virtual functions from SamplerBase. */
  //@{
  /**
   * Initialize the the sampler, possibly doing presampling of the
   * phase space.
   */
  virtual void initialize();

  /**
   * Generarate a new phase space point and return a weight associated
   * with it. This weight should preferably be 1.
   */
  virtual double generate();

  /**
   * Reject the last chosen phase space point.
   */
  virtual void rejectLast();

  /**
   * If the sampler is able to sample several different functions
   * separately, this function should return the last chosen
   * function. This default version always returns 0.
   */
  virtual int lastBin() const { return lastSampler ? lastSampler->bin() : 0; }

  /**
   * Return the total integrated cross section determined from the
   * Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSec() const {
    currentCrossSections();
    return theIntegratedXSec * nanobarn;
  }

  /**
   * Return the error on the total integrated cross section determined
   * from the Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSecErr() const {
    currentCrossSections();
    return theIntegratedXSecErr * nanobarn;
  }

  /**
   * Return the overestimated integrated cross section.
   */
  virtual CrossSection maxXSec() const {
    currentCrossSections();
    return abs(theIntegratedXSec) * nanobarn;
  }

  /**
   * Return the sum of the weights returned by generate() so far (of
   * the events that were not rejeted).
   */
  virtual double sumWeights() const {
    return theSumWeights;
  }
  //@}

protected:

  /**
   * Calculate cross sections from samplers, assuming the start of a
   * new iteration.
   */
  void updateCrossSections(bool firstTime = false);

  /**
   * Calculate cross sections from samplers at current state.
   */
  void currentCrossSections() const;

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


protected:

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();

private:

  /**
   * The bin sampler to use.
   */
  Ptr<BinSampler>::ptr theBinSampler;

  /**
   * Whether or not additional information should be printed to cout.
   */
  bool theVerbose;

  /**
   * True, if subprocesses should be selected flat. This is a debug
   * flag, cross section information and distributions will not be
   * correct.
   */
  bool theFlatSubprocesses;

  /**
   * True, if we are generating events.
   */
  bool isSampling;

  /**
   * The selector map for the bin samplers.
   */
  map<double,Ptr<BinSampler>::ptr> samplers;

  /**
   * The last selected bin sampler.
   */
  Ptr<BinSampler>::tptr lastSampler;

  /**
   * The number of events after which cross sections should truly be
   * updated. This is used to prevent exhaustive combination of
   * statistics when HepMC events are written out.
   */
  size_t theUpdateAfter;

  /**
   * The number of calls to currentCrossSections since the last
   * update.
   */
  mutable size_t crossSectionCalls;

  /**
   * True, if currentCrossSections has been called since the last call
   * to generate.
   */
  mutable bool gotCrossSections;

  /**
   * The integrated cross section in nanobarn
   */
  mutable double theIntegratedXSec;

  /**
   * The integrated cross section error in nanobarn
   */
  mutable double theIntegratedXSecErr;

  /**
   * The sum of weights
   */
  double theSumWeights;

  /**
   * The sum of absolute cross section.
   */ 
  double norm;

  /**
   * Map samplers to events to be skipped owing to encounter of a new
   * maximum.
   */
  map<Ptr<BinSampler>::tptr,unsigned long> skipMap;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralSampler & operator=(const GeneralSampler &);

};

}

#endif /* Herwig_GeneralSampler_H */
