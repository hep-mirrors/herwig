// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_Analysis2Base_H
#define Analysis2_Analysis2Base_H
//
// This is the declaration of the Analysis2Base class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Analysis2Base.fh"

#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"

#include "ThePEG/Utilities/StringUtils.h"

#include "Analysis2/Histogram2/Histogram2.h"
#include "Analysis2/Histogram2/Histogram2Output.h"
#include "EventExtractor.h"
#include "JetFinder.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 *
 * Analysis2 is the base class for all Analysis2
 * AnalysisHandlers. It uses the Histogram2 class
 * for histogramming and statistics and provides
 * helpers for setting up parallel runs.
 *
 * @see \ref Analysis2BaseInterfaces "The interfaces"
 * defined for Analysis2Base.
 */
class Analysis2Base: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Analysis2Base();

  /**
   * The destructor.
   */
  virtual ~Analysis2Base();
  //@}

public:

  /**
   * Insert an observable by building it from
   * a data file or a given bin range and bin number.
   * instring is taken to be a XML tag:
   * <options title="the title" xlabel="x label" ylabel="y label" \
   *  norm="..." datafile="file name" datatitle="name"/>
   * or
   * <options title="the title" xlabel="x label" ylabel="y label" \
   *  norm="..." lowerbound="..." upperbound="..." numbins="..."/>
   * Output options are set accordingly.
   */
  void insert (const string& name,
	       const string& instring,
	       Histogram2Options&);

  /**
   * Insert an observable.
   *
   * If not present, this will insert a channel named
   * Analysis2++ into the histogram. Events are booked to this
   * channel by default.
   */
  inline void insertObservable (const string& name,
				Histogram2Ptr,
				const Histogram2Options& =
				Histogram2Options());

  /**
   * Insert an observable by building it from a
   * data file. A channel Analysis2++ is inserted.
   */
  inline void insertObservable (const string& name,
				const string& dataName,
				const string& dataFile,
				const Histogram2Options& =
				Histogram2Options());

  /**
   * Return the histogram for the given
   * observable.
   */
  inline tcHistogram2Ptr histogram (const string&) const;

  /**
   * Book an event for the given observable.
   * The weight is looked up from the
   * event generator. The channel it is booked
   * to is named Analysis2++
   */
  void book (double evt, const string& name, double weight);

  /**
   * Normalization modes
   */
  enum NormalisationMode {
    FromMap = -1,
    NoNormalisation = 0,
    NormaliseToUnity = 1,
    NormaliseToXSec = 2,
    NormaliseToData = 3
  };

  /**
   * Finish an observable.
   *
   * If a channel Analysis2++ is found, this is normalized
   * depending on the normalization flag (see above).
   * If normalization to data is requested, but no data
   * channel is indicated, it is normalized to cross section.
   *
   * If we are counting per subprocess,
   * channels named Analysis2++-[Multiplicity] are summed up
   * to a channel named "Analysis2++", which is in turn normalized
   * as a single Analysis2++ channel.
   *
   * If a data channel name is indicated, the ratio
   * and chi2 are computed taking Analysis2++ as prediction.
   * The chi2/DOF, indexed by the observable's name, is written
   * to generator()->log().
   *
   */
  void finish (const string& name,
	       int normMode = FromMap,
	       bool combined = false);

  /**
   * Insert an observable by combining multiple runs.
   * @see Histogram2::combine
   */
  void combineObservable (const string& prefix,
			  const string& name,
			  unsigned int numRuns,
			  int normMode = FromMap);

  /**
   * A command indicating that we start combining
   * observables from parallel runs.
   */
  inline string startCombine (string);

  /**
   * A command to combine an observable from multiple
   * runs. Arguments need to be given as follows:
   *
   * [run prefix] [number of runs] [observable name] \
   * {normalization mode : none, unity, xsec, data}
   *
   * If observable = *, all present observables are
   * combined.
   */
  string combine (string);

  /**
   * A command indicating that we finished combining
   * observables from parallel runs.
   */
  inline string finishCombine (string);

  /**
   * Return the Histogram2Output object
   */
  inline tcHistogram2OutputPtr output () const;

  /**
   * Return the JetFinder
   */
  inline tcJetFinderPtr jetFinder () const;

  /**
   * Return the event extractor
   */
  inline tcEventExtractorPtr eventExtractor () const { return _eventExtractor; }

protected:

  /**
   * Access the histogram for the given
   * observable.
   */
  inline tHistogram2Ptr histogram (const string&);

  /**
   * Access the Histogram2Output object
   */
  inline tHistogram2OutputPtr output ();

  /**
   * Return the last presented event
   */
  inline tcEventPtr lastEvent () const;

  /**
   * Return the minimum over the presented
   * subprocess multiplicities.
   */
  inline unsigned int minMult () const;

  /**
   * Return the maximum over the presented
   * subprocess multiplicities.
   */
  inline unsigned int maxMult () const;

  /**
   * Access the JetFinder
   */
  inline tJetFinderPtr jetFinder ();

  /**
   * Access the event extractor
   */
  inline tEventExtractorPtr eventExtractor () { return _eventExtractor; }

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
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
  //@}

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

public:

  /**
   * Finalize the analysis. Called by dofinish.
   * Dumps the histograms to h2 files for later recovery.
   */
  inline void finalize ();

protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();

  //@}

private:

  /**
   * The output used for histograms
   */
  Histogram2OutputPtr _output;

  /**
   * The observables mapped by their name.
   */
  map<string,Histogram2Ptr> _histograms;

  /**
   * The observables mapped to their output
   * options.
   */
  map<string,Histogram2Options> _outputOptions;

  /**
   * Observables mapped to data channels
   */
  map<string,string> _datachannels;

  /**
   * Observables mapped to normalisation mode
   */
  map<string,int> _normalisation;

  /**
   * Wether we should book a channel per
   * hard subprocess multiplitiplicity.
   */
  bool _bookPerSubprocess;

  /**
   * The last presented events.
   */
  tcEventPtr _lastEvent;

  /**
   * The minimum multiplicity of the
   * hard processes presented.
   */
  unsigned int _minMult;

  /**
   * The maximum multiplicity of the
   * hard processes presented.
   */
  unsigned int _maxMult;

  /**
   * Wether or not we're part of a parallel run.
   */
  bool _parallel;

  /**
   * The jet finder to be used
   */
  JetFinderPtr _jetFinder;

  /**
   * The event extractor to be used
   */
  EventExtractorPtr _eventExtractor;

  /**
   * Backup parallel flag
   */
  bool _backParallel;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<Analysis2Base> initAnalysis2Base;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Analysis2Base & operator=(const Analysis2Base &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Analysis2Base. */
template <>
struct BaseClassTrait<Analysis2::Analysis2Base,1> {
  /** Typedef of the first base class of Analysis2Base. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Analysis2Base class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::Analysis2Base>
  : public ClassTraitsBase<Analysis2::Analysis2Base> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::Analysis2Base"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Analysis2Base is implemented. It may also include several, space-separated,
   * libraries if the class Analysis2Base depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#include "Analysis2Base.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Analysis2Base.tcc"
#endif

#endif /* Analysis2_Analysis2Base_H */
