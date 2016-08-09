// -*- C++ -*-
//
// GeneralSampler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
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
  virtual int lastBin() const { return lastSampler() ? lastSampler()->bin() : 0; }

  /**
   * Return the total integrated cross section determined from the
   * Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSec() const {
    currentCrossSections();
    return theIntegratedXSec;
  }

  /**
   * Return the error on the total integrated cross section determined
   * from the Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSecErr() const {
    currentCrossSections();
    return theIntegratedXSecErr;
  }

  /**
   * Return the overestimated integrated cross section.
   */
  virtual CrossSection maxXSec() const { 
    if ( theAddUpSamplers )
      return SamplerBase::maxXSec();
    return theMaxWeight*nanobarn;
  }

  /**
   * Return the sum of the weights returned by generate() so far (of
   * the events that were not rejeted).
   */
  virtual double sumWeights() const { return theSumWeights; }

  /**
   * Return the sum of the weights squaredreturned by generate() so far (of
   * the events that were not rejeted).
   */
  virtual double sumWeights2() const { return theSumWeights2; }

  /**
   * Return the number of attempts
   */
  virtual double attempts() const { 
    if ( theAddUpSamplers )
      return SamplerBase::attempts();
    return theAttempts;
  }

  /**
   * Return the number of accepts
   */
  double accepts() const { return theAccepts; }

  //@}

  /**
   * Return the samplers
   */
  const map<double,Ptr<BinSampler>::ptr>& samplers() const { return theSamplers; }

  /**
   * Return the bin sampler
   */
  Ptr<BinSampler>::ptr binSampler() const { return theBinSampler; }

  /**
   * Return the last selected bin sampler
   */
  Ptr<BinSampler>::tptr lastSampler() const { return theLastSampler; }

  /**
   * True if we should do weighted events
   */
  bool weighted() const { return eventHandler()->weighted(); }


  /** 
   * True if the sampler runs in Allmostunweighted mode.
   */ 

  bool almostUnweighted() const { return theAlmostUnweighted; }

public:

  /**
   * Return the XML element containing the grids
   */
  const XML::Element& grids() const { return theGrids; }

  /**
   * Access the XML element containing the grids
   */
  XML::Element& grids() { return theGrids; }

  /**
   * Write out grids
   */
  void writeGrids() const;

  /**
   * Read in grids
   */
  void readGrids();
  
  /**
   * Return the number of integration jobs which were actually created.
   */
  unsigned int integrationJobsCreated() {
    return theIntegrationJobsCreated;
  }

  /**
   * An external hook to prepare the sampler for generating events, e.g. by
   * combining grid files from parallel integration runs.
   */
  virtual void prepare();

protected:

  /**
   * Access the samplers
   */
  map<double,Ptr<BinSampler>::ptr>& samplers() { return theSamplers; }

  /**
   * Set the last selected bin sampler
   */
  void lastSampler(Ptr<BinSampler>::tptr s) { theLastSampler = s; }

  /**
   * Calculate cross sections from samplers at current state.
   */
  void currentCrossSections() const;

  /**
   * Update the sampler selection
   */
  void updateSamplers();

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

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
   * Whether or not additional information should be printed to cout.
   */
  bool theVerbose;

  /**
   * The XML element containing the grids
   */
  XML::Element theGrids;

  /**
   * The bin sampler to use.
   */
  Ptr<BinSampler>::ptr theBinSampler;

  /**
   * The selector map for the bin samplers.
   */
  map<double,Ptr<BinSampler>::ptr> theSamplers;

  /**
   * The last selected bin sampler.
   */
  Ptr<BinSampler>::tptr theLastSampler;

  /**
   * The integrated cross section
   */
  mutable CrossSection theIntegratedXSec;

  /**
   * The integrated cross section error
   */
  mutable CrossSection theIntegratedXSecErr;

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
   * The sum of weights
   */
  double theSumWeights;

  /**
   * The sum of weights squared
   */
  double theSumWeights2;

  /**
   * The number of attempts
   */
  double theAttempts;

  /**
   * The number of accepts
   */
  double theAccepts;

  /**
   * The maximum weight encountered
   */
  double theMaxWeight;

  /**
   * True, if cross sections are to be combined from each sampler
   * individually
   */
  bool theAddUpSamplers;

  /**
   * True, if the global maximum weight should be used as
   * reference. If not, the maximum weights of individual samplers are
   * used, and selection probabilities fro the samplers are adjusted
   * accordingly.
   */
  bool theGlobalMaximumWeight;

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
   * A minimum selection probability for each sampler
   */
  double theMinSelection;

  /**
   * True, if information for combining unnormalized runs should be
   * printed out
   */
  bool runCombinationData;

  /**
   * True, if we should perform an almost unweighted sampling
   */
  bool theAlmostUnweighted;

  /**
   * Number of points which exceeded the maximum
   */
  unsigned long maximumExceeds;

  /**
   * The average relative deviation from the maximum weight
   */
  double maximumExceededBy;

  /**
   * The correct cross section as one would exspect with
   * almostUnweighted. 
   */

  double correctWeights;

  /**
   * Enhancement factor to the maximum weight.
   * This is to get less maximumExceeds. 
   */

  double theMaxEnhancement;

  /**
   * True, if grids have already been read.
   */
  bool didReadGrids;

  /**
   * True, if parallel subprocess integration should be enabled
   */
  bool theParallelIntegration;

  /**
   * The number of subprocesses to integrate per job
   */
  unsigned int theIntegratePerJob;

  /**
   * The maximum number of integration jobs to be created
   */
  unsigned int theIntegrationJobs;

  /**
   * The number of integration jobs which were actually created
   */
  unsigned int theIntegrationJobsCreated;
  
  /**
   * Indicate that initialization is only reading a grid.
   */
  bool justAfterIntegrate;

  /**
   * True, if grids should be written at the end of a run
   */
  bool theWriteGridsOnFinish;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralSampler & operator=(const GeneralSampler &);

};

}

#endif /* Herwig_GeneralSampler_H */
