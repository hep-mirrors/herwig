// -*- C++ -*-
#ifndef Herwig_HardProcessAnalysis_H
#define Herwig_HardProcessAnalysis_H
//
// This is the declaration of the HardProcessAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HardProcessAnalysis class.
 *
 * @see \ref HardProcessAnalysisInterfaces "The interfaces"
 * defined for HardProcessAnalysis.
 */
class HardProcessAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HardProcessAnalysis();

  /**
   * The destructor.
   */
  virtual ~HardProcessAnalysis();
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardProcessAnalysis & operator=(const HardProcessAnalysis &) = delete;

  /**
   * Differential information per outgoing parton
   */
  struct Histograms {

    /**
     * The constructor
     */
    Histograms() {}

    /**
     * The constructor
     */
    explicit Histograms(Energy ECM, unsigned int theNBins);

    /**
     * Analyse given momentum
     */
    void fill(const Lorentz5Momentum& p, double weight);

    /**
     * Finalize given process id and cross section.
     */
    void finalize(ostream& dat,
		  ostream& plot,
		  const string& subpro,
		  size_t legid,
		  double norm,
		  bool theUnitWeights);

    /**
     * Pt spectrum
     */
    HistogramPtr transverse;

    /**
     * Rapidity distribution
     */
    HistogramPtr rapidity;

    /**
     * Azimuthal angle distribution
     */
    HistogramPtr phi;

  };

  /**
   * Outgoing partons and x distributions
   */
  struct AllHistograms {

    /**
     * Outgoing partons
     */
    vector<Histograms> outgoing;

    /**
     * x1 distribution
     */
    HistogramPtr x1;

    /**
     * x2 distribution
     */
    HistogramPtr x2;

    /**
     * sqrt(shat) distribution
     */
    HistogramPtr sshat;

    /**
     * y distribution
     */
    HistogramPtr rapidity;

    /**
     * The sum of weights
     */
    double sumWeights;

  };

  /**
   * Histograms per subprocess
   */
  map<vector<string>,AllHistograms> histogramData;

  /**
   * The total sum of weights
   */
  double sumWeights;

  /**
   * Analyze a given final state
   */
  void fill(PPair, ParticleVector, double);

  /**
   * The number of bins to use
   */
  unsigned int theNBins;

  /**
   * True, if unit weights should be booked
   */
  bool theUnitWeights;

  /**
   * True, if subprocesses should be distinguished by initial state
   */
  bool theSplitInitialStates;

  /**
   * True, if partons should be handled as jets irrespective of flavour
   */
  bool thePartonsAreJets;

};

}

#endif /* Herwig_HardProcessAnalysis_H */
