// -*- C++ -*-
#ifndef HERWIG_LPairAnalysis_H
#define HERWIG_LPairAnalysis_H
//
// This is the declaration of the LPairAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The LPairAnalysis class is designed to analyse lepton pairs.
 * Therefore the lepton and antilepton with highest pt are selected from
 * the event and some observables are computed and booked into
 * histograms.  The analysis is intended to anaylse top pair events
 * where both W's are decaying semileptonically.  Therefore, only the
 * semileptonic top quark decay channels should be switched on.  If no
 * pair of (opposite sign) leptons is found in the final state a warning
 * message is printed and nothing is booked.  Output is created as
 * topdrawer file
 *
 * @see \ref LPairAnalysisInterfaces "The interfaces"
 * defined for LPairAnalysis.
 */
class LPairAnalysis: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  LPairAnalysis();

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

public:

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class without persistent data.
   */
  static NoPIOClassDescription<LPairAnalysis> initLPairAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LPairAnalysis & operator=(const LPairAnalysis &);

private:
  /**
   *   \f$p_T\f$ of the leptons
   */
  Histogram _ptp;
  Histogram _ptm;
  Histogram _ptpair;

  /**
   *   \f$E_T\f$ of the leptons
   */
  Histogram _etp;
  Histogram _etm;
  Histogram _etpair;

  /**
   *   Energy of the leptons
   */
  Histogram _ep;
  Histogram _em;
  Histogram _epair;

  /**
   *  Rapidity of the leptons
   */
  Histogram _rapp;
  Histogram _rapm;
  Histogram _rappair;

  /**
   *  Azimuth of the leptons
   */
  Histogram _phip;
  Histogram _phim;
  Histogram _deltaphi;
  
  /**
   *  Invariant mass of the pair
   */
  Histogram _mpair;

  /**
   *  scalar sums of Et, pt
   */
  Histogram _etsum;
  Histogram _ptsum;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LPairAnalysis. */
template <>
struct BaseClassTrait<Herwig::LPairAnalysis,1> {
  /** Typedef of the first base class of LPairAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LPairAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LPairAnalysis>
  : public ClassTraitsBase<Herwig::LPairAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LPairAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the LPairAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LPairAnalysis_H */
