// -*- C++ -*-
//
// BinSampler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_BinSampler_H
#define Herwig_BinSampler_H
//
// This is the declaration of the BinSampler class.
//

#include "ThePEG/Handlers/StandardEventHandler.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Repository/UseRandom.h"

#include "MultiIterationStatistics.h"
#include "Remapper.h"

namespace Herwig {

using namespace ThePEG;

class GeneralSampler;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief BinSampler samples XCombs bins. This default implementation
 * performs flat MC integration.
 *
 * @see \ref BinSamplerInterfaces "The interfaces"
 * defined for BinSampler.
 */
class BinSampler: public Herwig::MultiIterationStatistics {

public:

  /**
   * The default constructor.
   */
  BinSampler();

public:

  /**
   * Clone this object.
   */
  Ptr<BinSampler>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<BinSampler>::ptr>(clone());
  }

public:

  /**
   * Evaluate the cross section
   */
  double evaluate(vector<double> p,
		  bool remap = true);

  /**
   * Return the bias with which this sampler is selected. The sampler
   * needs to divide out this bias in its weight calculation.
   */
  double bias() const { return theBias; }

  /**
   * Set the bias with which this sampler is selected.
   */
  void bias(double b) { theBias = b; }

  /**
   * Set the event handler
   */
  void eventHandler(tStdEHPtr eh) { theEventHandler = eh; }

  /**
   * Return the event handler
   */
  tStdEHPtr eventHandler() const { return theEventHandler; }

  /**
   * Set the containing sampler
   */
  void sampler(Ptr<GeneralSampler>::tptr);

  /**
   * Get the containing sampler
   */
  Ptr<GeneralSampler>::tptr sampler() const;

  /**
   * Return the bin
   */
  int bin() const { return theBin; }

  /**
   * Set the bin
   */
  void bin(int b) { theBin = b; }

  /**
   * Return a string describing the process handled by this sampler.
   */
  string process() const;

  /**
   * Return a short string describing the process handled by this sampler.
   */
  string shortprocess() const;
  
  /**
   * Return a string identifying the process handled by this sampler.
   */
  string id() const;

  /**
   * Return the last generated point.
   */
  const vector<double>& lastPoint() const { return theLastPoint; }

  /**
   * Access the last generated point.
   */
  vector<double>& lastPoint() { return theLastPoint; }

  /**
   * Return the reference weight to be used
   */
  double referenceWeight() const { return theReferenceWeight; }

  /**
   * Set the reference weight to be used
   */
  void referenceWeight(double w) { theReferenceWeight = w; }

  /**
   * Return true, if this sampler can provide unweighted events; if
   * the proposal density is not an overestimate, weights larger than
   * one can be generated, the handling of these points being subject
   * to the GeneralSampler class.
   */
  virtual bool canUnweight() const { return true; }

  /**
   * Return true, if this sampler adapts on the fly while generating
   * events. Cross sections in the GeneralSampler class are calculated
   * from adding up the cross sections quoted by individual samplers.
   */
  virtual bool adaptsOnTheFly() const { return false; }

  /**
   * If this sampler features a compensation algorithm, return true if
   * more events need to be generated to finish the compensation.
   */
  virtual bool compensating() const { return false; }

  /**
   * Return true, if weighted events should be generated
   */
  bool weighted() const { return theWeighted; }

  /**
   * Indicate that weighted events should be generated
   */
  void doWeighted(bool yes = true) { theWeighted = yes; }

  /**
   * Exception to be thrown if cross section information should be updated.
   */
  struct NextIteration {};

  /**
   * Generate the next point and return its weight; store the point in
   * lastPoint().
   */
  virtual double generate();

  /**
   * Fill and finalize the remappers present
   */
  void fillRemappers(bool progress);

  /**
   * Write remappers to grid file
   */
  void saveRemappers() const;

  /**
   * Write integration data to grid files
   */
  void saveIntegrationData() const;

  /**
   * Save grid data
   */
  virtual void saveGrid() const {}

  /**
   * Read integration data from grid files
   */
  void readIntegrationData();

  /**
   * Read remappers from grid file
   */
  void setupRemappers(bool progress);

  /**
   * Run a single iteration of n points, optionally printing a
   * progress bar to cout. Calls generate n times.
   */
  void runIteration(unsigned long n, bool progress);

  /**
   * Adapt this sampler after an iteration has been run
   */
  virtual void adapt() {}

  /**
   * Initialize this bin sampler. This default version calls runIteration.
   */
  virtual void initialize(bool progress);

  /**
   * Return true, if this sampler has already been initialized.
   */
  bool initialized() const { return theInitialized; }

  /**
   * Indicate that this sampler has already been initialized.
   */
  void isInitialized() { theInitialized = true; }

  /**
   * Return true, if integration has already been performed
   */
  bool integrated() const { return theIntegrated; }

  /**
   * Return true, if remappers have been set up
   */
  bool remappersFilled() const { return theRemappersFilled; }

  /**
   * Return true, if grid data exists for this sampler.
   */
  virtual bool existsGrid() const { return false; }

  /**
   * Return true, if this sampler has already read grid data.
   */
  bool hasGrids() const { return theHasGrids; }

  /**
   * Indicate that this sampler has already read grid data.
   */
  void didReadGrids() { theHasGrids = true; }

  /**
   * Finalize this sampler.
   */
  virtual void finalize(bool); 

  /**
   * Return the total integrated cross section determined from the
   * Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSec() const {
    return averageWeight()*nanobarn;
  }

  /**
   * Return the error on the total integrated cross section determined
   * from the Monte Carlo sampling so far.
   */
  virtual CrossSection integratedXSecErr() const {
    return sqrt(abs(averageWeightVariance()))*nanobarn;
  }

   /**
   * Define the key for the collinear subtraction data.
   */

  
  struct RandomNumberHistogram {

    /**
     * The lower bound
     */
    double lower;

    /**
     * The bins, indexed by upper bound.
     */
    map<double,double > bins;
    
    map<double,double > binsw1;
    /**
     * Constructor
     */
    RandomNumberHistogram(double low = 0.0, 
			 double up = 1., 
			 unsigned int nbins = 20);

    /**
     * Book an event.
     */
    void book(double inv, double weight) {
      map<double,double>::iterator b =	bins.upper_bound(inv);
      if ( b == bins.end() ) return;
      b->second = b->second+weight;
      map<double,double>::iterator b2 =	binsw1.upper_bound(inv);
      if ( b2 == binsw1.end() ) return;
      b2->second = b2->second+1.;

    }

    /**
     * Write to file given name and invariant.
     */
    void dump(const std::string& folder,const std::string& prefix, const std::string& process,const int NR)const;


  };
     
  typedef pair<string,size_t > RandomNumberIndex;
  
  map<RandomNumberIndex,pair<RandomNumberHistogram,double> > RandomNumberHistograms;

public:

  /**
   * Return the dimension.
   */
  int dimension() const { return theEventHandler->nDim(bin()); }

  /**
   * Return the number of points to be used for initial integration.
   */
  unsigned long initialPoints() const { return theInitialPoints; }

  /**
   * Set the number of points to be used for initial integration.
   */
  void initialPoints(unsigned long n) { theInitialPoints = n; }

  /**
   * Return the number of iterations to be considered for initialization.
   */
  size_t nIterations() const { return theNIterations; }

  /**
   * Set the number of iterations to be considered for initialization.
   */
  void nIterations(size_t n) { theNIterations = n; }

  /**
   * Set the factor to enhance the number of points for the next
   * iteration.
   */
  void enhancementFactor(double f) { theEnhancementFactor = f; }

  /**
   * Return the factor to enhance the number of points for the next
   * iteration.
   */
  double enhancementFactor() const { return theEnhancementFactor; }
  
  /**
   * Return the folder for the random number plots.
   */
  string randomNumberString() const {return theRandomNumbers;}

  /**
   * In the AlmostUnweighted mode we do not need to unweight 
   * the events to the reference weight. 
   * Kappa reduces effectivly the reference weight.
   * This can be used for processes, where unweighting 
   * is hardly feasable.
   */
  double kappa() const {return theKappa;}

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
   * The bias with which this sampler is selected.
   */
  double theBias;

  /**
   * True, if weighted events should be generated
   */
  bool theWeighted;

  /**
   * The number of points to use for initial integration.
   */
  unsigned long theInitialPoints;

  /**
   * The number of iterations to be considered for initialization.
   */
  size_t theNIterations;

  /**
   * Factor to enhance the number of points for the next iteration.
   */
  double theEnhancementFactor;

  /**
   * Switch to count only non zero weights in presampling.
   */

  bool theNonZeroInPresampling; 
  
  /**
   * Switch to require that we get half of the points 
   * in each iteration below the maximum weight of the iteration.
   */    
    
  bool theHalfPoints;

  /**
   * The maximum number of allowed new maxima, 
   * in combination with HalfPoints, in order to prevent unstable
   * processes.
   */
  int theMaxNewMax;

  /**
   * The reference weight to be used
   */
  double theReferenceWeight;

  /**
   * The bin to be sampled.
   */
  int theBin;

  /**
   * Wether or not this sampler has already been initialized.
   */
  bool theInitialized;

  /**
   * The last generated point.
   */
  vector<double> theLastPoint;

  /**
   * The event handler to be used.
   */
  tStdEHPtr theEventHandler;

  /**
   * The containing sampler
   */
  Ptr<GeneralSampler>::tptr theSampler;
  
  /**
   * Folder for the random number plots.
   */
  string theRandomNumbers;

  /**
   * Remapper objects indexed by dimension
   */
  map<size_t,Remapper> remappers;

  /**
   * The number of points to be used for initial filling of the remappers
   */
  unsigned long theRemapperPoints;

  /**
   * True if channels should get a remapper
   */
  bool theRemapChannelDimension;

  /**
   * The number of bins to be used for luminosity dimensions
   */
  unsigned long theLuminosityMapperBins;

  /**
   * The number of bins to be used for any other dimension
   */
  unsigned long theGeneralMapperBins;

  /**
   * The minimum selection probability for remapper bins
   */
  double theRemapperMinSelection;

  /**
   * True, if integration has already be performed
   */
  bool theIntegrated;

  /**
   * True, if remappers have been set up
   */
  bool theRemappersFilled;

  /**
   * True, if this sampler has already read grid data.
   */
  bool theHasGrids;



  /**
   * In the AlmostUnweighted mode we do not need to unweight 
   * the events to the reference weight. 
   * Kappa reduces effectivly the reference weight.
   * This can be used for processes, where unweighting 
   * is hardly feasable.
   */
  double theKappa;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BinSampler & operator=(const BinSampler &) = delete;

};

}

#endif /* Herwig_BinSampler_H */
