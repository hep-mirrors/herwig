// -*- C++ -*-
#ifndef HERWIG_DrellYanPT_H
#define HERWIG_DrellYanPT_H
//
// This is the declaration of the DrellYanPT class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * 
 * This AnalysisHandler books histograms
 * of the weak boson's pt in Drell-Yan
 * processes.
 *
 * An Interface switch allows you to choose
 * between output of the histograms
 * as TopDraw or plain ASCII files, which
 * may be processed further, e.g. by gnuplot.
 *
 * @see \ref DrellYanPTInterfaces "The interfaces"
 * defined for DrellYanPT.
 */
class DrellYanPT: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  DrellYanPT();

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
  //@}

private:

  /**
   * Histogram of the Z's pt
   */
  Histogram _Zpt;

  /**
   * Histogram of the W+'s pt
   */
  Histogram _Wppt;

  /**
   * Histogram of the W-'s pt
   */
  Histogram _Wmpt;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DrellYanPT & operator=(const DrellYanPT &) = delete;

};

}

#endif /* HERWIG_DrellYanPT_H */
