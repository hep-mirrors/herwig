// -*- C++ -*-
#ifndef HERWIG_TauTo3MesonAnalysis_H
#define HERWIG_TauTo3MesonAnalysis_H
//
// This is the declaration of the TauTo3MesonAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The TauTo3MesonAnalysis class is designed to perform the analysis of the
 *  mass distribution of the hadronic decay products of the \f$\tau\f$ in the decays
 * - \f$\tau^-\to\nu_tau \pi^+\pi^-\pi^-  \f$
 * - \f$\tau^-\to\nu_tau \pi^0\pi^0\pi^-  \f$
 * - \f$\tau^-\to\nu_tau K^-K^+\pi^-      \f$
 * - \f$\tau^-\to\nu_tau K^0\bar{K}^0\pi^-\f$
 * - \f$\tau^-\to\nu_tau K^-K^0\pi^0      \f$
 * - \f$\tau^-\to\nu_tau \pi^0-\pi^0K^-   \f$
 * - \f$\tau^-\to\nu_tau K^-\pi^-\pi^+    \f$
 * - \f$\tau^-\to\nu_tau \pi^-K^0\pi^0    \f$
 * - \f$\tau^-\to\nu_tau \pi^-\pi^0\eta   \f$
 * - \f$\tau^-\to\nu_tau \pi^-\pi^0\gamma \f$
 *
 * @see \ref TauTo3MesonAnalysisInterfaces "The interfaces"
 * defined for TauTo3MesonAnalysis.
 */
class TauTo3MesonAnalysis: public AnalysisHandler {

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

  using AnalysisHandler::analyze;
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

protected:

  /** @name Standard Interfaced functions. */
  //@{
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
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<TauTo3MesonAnalysis> initTauTo3MesonAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TauTo3MesonAnalysis & operator=(const TauTo3MesonAnalysis &);

private:

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^+\pi^-\pi^-  \f$
   */
  vector<HistogramPtr> _m3pippimpim;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^0\pi^0\pi^-  \f$
   */
  vector<HistogramPtr> _m3pi0pi0pim;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^-K^+\pi^-      \f$
   */
  vector<HistogramPtr> _m3kmpimkp;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^0\bar{K}^0\pi^-\f$
   */
  vector<HistogramPtr> _m3k0pimk0;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^-K^0\pi^0      \f$
   */
  vector<HistogramPtr> _m3kmpi0k0;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^0\pi^0K^-   \f$
   */
  vector<HistogramPtr> _m3pi0pi0km;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^-\pi^-\pi^+    \f$
   */
  vector<HistogramPtr> _m3kmpimpip; 

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^-K^0\pi^0    \f$
   */
  vector<HistogramPtr> _m3pimk0pi0;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^-\pi^0\eta   \f$
   */
  vector<HistogramPtr> _m3pimpi0eta;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau \pi^-\pi^0\gamma \f$
   */
  vector<HistogramPtr> _m3pimpi0gamma;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^0_SK^0_S\pi^-\f$
   */
  vector<HistogramPtr> _m3kspimks;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^0_LK^0_L\pi^-\f$
   */
  vector<HistogramPtr> _m3klpimkl;

  /**
   *  Histograms for \f$\tau^-\to\nu_tau K^0_SK^0_L\pi^-\f$
   */
  vector<HistogramPtr> _m3kspimkl;


};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TauTo3MesonAnalysis. */
template <>
struct BaseClassTrait<Herwig::TauTo3MesonAnalysis,1> {
  /** Typedef of the first base class of TauTo3MesonAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TauTo3MesonAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TauTo3MesonAnalysis>
  : public ClassTraitsBase<Herwig::TauTo3MesonAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TauTo3MesonAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TauTo3MesonAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class TauTo3MesonAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwTauAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_TauTo3MesonAnalysis_H */
