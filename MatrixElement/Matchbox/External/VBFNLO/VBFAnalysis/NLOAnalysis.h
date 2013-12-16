// -*- C++ -*-
#ifndef ThePEG_NLOAnalysis_H
#define ThePEG_NLOAnalysis_H
//
// This is the declaration of the NLOAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "BasicHistogramCollection.h"
#include "arHistogram.fh"
#include "arHistogram.h"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the NLOAnalysis class.
 *
 * @see \ref NLOAnalysisInterfaces "The interfaces"
 * defined for NLOAnalysis.
 */
class NLOAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLOAnalysis();

  /**
   * The destructor.
   */
  virtual ~NLOAnalysis();
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

  /**
   * Return a LorentzTransform which would put the event in the
   * desired Lorentz frame.
   * @param event a pointer to the Event to be considered.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tcEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the weight of the current event.
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight the weight of the current event.
   */
  virtual void analyze(tPPtr particle, double weight);
  //@}

  /**
   * This implements the actual analysis and is called with different
   * sets of particles and weights depending if the event is real emission
   * or subtraction or born-type
   */
  void doAnalyze(tEventPtr event, const vector<fastjet::PseudoJet>& input_particles,
		 const vector<fastjet::PseudoJet>& hard_partons, const double w);

  /**
   * Check if the particle is allowed for jet reconstruction
   */
  bool allowedParticle (const Particle &p);

  /**
   * Check for cuts
   */
  bool passCuts(vector<fastjet::PseudoJet>);

  /**
   * Push particles into a PseudoJet vector with the option of sorting out
   */
  vector<fastjet::PseudoJet> recombinables(const ParticleVector& p, bool sort_out = true);

  /**
   * Push particles into a PseudoJet vector with the option of sorting out
   */
  vector<fastjet::PseudoJet> recombinables(const tPVector& p, bool sort_out = true);

  /**
   * Recombine to jets using fastjet methods
   */
  vector<fastjet::PseudoJet> recombine(const vector<fastjet::PseudoJet>& p);

  /**
   * Return jets within the detector range
   */
  vector<fastjet::PseudoJet> getInRange(const vector<fastjet::PseudoJet>& j);

  /**
   * Get the Higgs
   */
  fastjet::PseudoJet getHiggs(const ParticleVector& p);

  /**
   * Set the cached jets
   */
  void setJetCache(const vector<fastjet::PseudoJet>& j) {
    theJetCacheSet = true;
    jetCache = j; 
  }

  /**
   * Clear the cached jets
   */
  void clearJetCache() { 
    theJetCacheSet = false;
    jetCache = vector<fastjet::PseudoJet>(); 
  }

  /**
   * Flag to see if the cache is set
   */
  bool cacheIsSet() { return theJetCacheSet; }

  /**
   * Flag to see if running PlainNLOmode
   */
  bool runsPlainNLO() { return plainNLO; }

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

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   *
   * This method is called just before
   * running starts through the initrun() method to indicate that the
   * actual running is to start. When implemented by a sub class it is
   * important that the doinitrun() method of the base class is called
   * first and then, if the initialization of this object depends on
   * other objects, that the initrun() method of these objects are
   * called. Only then should the class-local initialization
   * proceed. To avoid circular loops, it is important that the
   * doinitrun() method is called for the base class, while the
   * initrun() method is called for other objects.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   *
   * This method is called after the running
   * phase through the finish() and can eg. be used to write out
   * statistics. When implemented by a sub class it is important that
   * the dofinish() method of the base class is called while the
   * finish() methd is called for other objects.
   */
  virtual void dofinish();


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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NLOAnalysis & operator=(const NLOAnalysis &);

  /**
   * The histograms of the events
   */
  arnold::BasicHistogramCollection events;

  /**
   * The histograms for the higgs
   */
  arnold::SingleJetHistogramCollection higgs;

  /**
   * The detector range in pseudorapidity
   */
  double theEtaDetector;
  
  /**
   * The jet algorithm cone parameter
   */
  double theRparam;

  /**
   * The jet algoritm minimum pt definition
   */
  double theJetDefPt;

  /**
   * Set to true for a plain NLO calculation
   */
  bool plainNLO;

  /**
   * The cached jets
   */
  vector<fastjet::PseudoJet> jetCache;

  /**
   * The jets
   */
  vector<fastjet::PseudoJet> jets;

  /**
   * Flag to check if the jet cache is filled
   */
  bool theJetCacheSet;

  bool debuginfo;

  double theRjj_min;
  double theRapidityGap;
  double themjj_min;
  bool theopposite_dir;

  unsigned int thenjets_min, thenjets_ex;
  double theDelYhj_min;

};

}

#endif /* ThePEG_NLOAnalysis_H */
