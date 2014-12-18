// -*- C++ -*-
#ifndef Herwig_WeightAnalyzer_H
#define Herwig_WeightAnalyzer_H
//
// This is the declaration of the WeightAnalyzer class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the WeightAnalyzer class.
 *
 * @see \ref WeightAnalyzerInterfaces "The interfaces"
 * defined for WeightAnalyzer.
 */
class WeightAnalyzer: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  WeightAnalyzer();

  /**
   * The destructor.
   */
  virtual ~WeightAnalyzer();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //@}

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
   * The sum of weights
   */
  double sumWeights;

  /**
   * The sum of positive weights
   */
  double sumPositiveWeights;

  /**
   * The sum of negative weights
   */
  double sumNegativeWeights;

  /**
   * The sum of weights calculated by subprocess group weights
   */
  double sumGroupWeights;

  /**
   * The sum of positive weights calculated by subprocess group weights
   */
  double sumPositiveGroupWeights;

  /**
   * The sum of negative weights calculated by subprocess group weights
   */
  double sumNegativeGroupWeights;

  /**
   * The maximum deviation of the group weight sum from one
   */
  double maxDeviationGroupWeight;

  /**
   * The maximum deviation of the event weight sum from the overall
   * event weight
   */
  double maxDeviationEventWeight;

  /**
   * Total number of positive weights
   */
  double nPositiveWeights;

  /**
   * Total number of negative weights
   */
  double nNegativeWeights;

  /**
   * The maximum postive weight
   */
  double maxPositiveWeight;

  /**
   * The maximum absolute negative weight
   */
  double maxNegativeWeight;

  /**
   * Histogram of positive weight occurences
   */
  map<double,double> positiveWeightDistribution;

  /**
   * Histogram of negative weight occurences
   */
  map<double,double> negativeWeightDistribution;
  
  
  /**
   * Gnuplot output
   */
  bool gnuplot;
  
  

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WeightAnalyzer & operator=(const WeightAnalyzer &);

};

}

#endif /* Herwig_WeightAnalyzer_H */
