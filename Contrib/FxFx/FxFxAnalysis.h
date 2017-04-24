// -*- C++ -*-
#ifndef THEPEG_FxFxAnalysis_H
#define THEPEG_FxFxAnalysis_H
//
// This is the declaration of the FxFxAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Rivet/AnalysisHandler.hh"
#include "FxFxReader.h"
#include "FxFxEventHandler.h"

namespace ThePEG {

/**
 * Here is the documentation of the FxFxAnalysis class.
 *
 * @see \ref FxFxAnalysisInterfaces "The interfaces"
 * defined for FxFxAnalysis.
 */
class FxFxAnalysis: public ThePEG::AnalysisHandler {

public:
  
  /**
   * The default constructor.
   */
  FxFxAnalysis();
  
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
  virtual void analyze(ThePEG::tEventPtr event, long ieve, int loop, int state);

   /**
   * Produca a HepMC event for the given subprocess
   */
  HepMC::GenEvent * makeEvent(tEventPtr event, tSubProPtr sub, long no,
			      Energy eUnit, Length lUnit,
			      CrossSection xsec, CrossSection xsecErr) const;

     /**
   * Produca a HepMC event for the given subprocess
   */
  HepMC::GenEvent * makeEventW(tEventPtr event, tSubProPtr sub, long no,
					  Energy eUnit, Length lUnit, 
                                           CrossSection xsec, CrossSection xsecErr, double evoptweight, double centralweight) const;

  int _i;

  int _numweights;
    
  //@}
  
public:
  
  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(ThePEG::PersistentOStream & os) const;
  
  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(ThePEG::PersistentIStream & is, int version);
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
  virtual ThePEG::IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual ThePEG::IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the read phase.
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
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ThePEG::ClassDescription<FxFxAnalysis> initFxFxAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FxFxAnalysis & operator=(const FxFxAnalysis &);

private:
  /**
   * The PDG ID to be used for remnants
   */
  long _remnantId;

  /**
   *  The HepMC format
   */
  int _format;

  /**
   * Selector for the choice of units
   */
  int _unitchoice;

  /**
   * Choice of output precision in GenEvent format
   */
  unsigned int _geneventPrecision;

  /**
   *  The Analyses to use
   */
  vector<string> _analyses;

  /**
   * The base name of the output file.
   */
  string filename;

  /**
   * Enable debugging information from FxFx
   */
  bool debug;

  /** 
   * Enable use of optional weights in analysis
   */  
  bool useoptweights = false;

   /** 
   * normalize optional weights to the central weight
   */  

  bool normoptweights = false;
  
  /**
   *  The FxFxAnalysisHandler
   */
  Rivet::AnalysisHandler * _rivet;

     /**
   *  The FxFxAnalysisHandlers for multiple weights
   */
  Rivet::AnalysisHandler * _rivetMULTI[120];


  /** 
   * holders of weights and cross section
   */

  std::vector< std::pair<int,double> > OptWeights;
  std::vector<double> OptXS;
  map<string,CrossSection> optxsec;


  /**
   *  Event count
   */
  unsigned long _nevent;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FxFxAnalysis. */
template <>
struct BaseClassTrait<FxFxAnalysis,1> {
  /** Typedef of the first base class of FxFxAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FxFxAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<FxFxAnalysis>
  : public ClassTraitsBase<FxFxAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::FxFxAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FxFxAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class FxFxAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "FxFxAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_FxFxAnalysis_H */
