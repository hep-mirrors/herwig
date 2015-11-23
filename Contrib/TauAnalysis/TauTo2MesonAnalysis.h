// -*- C++ -*-
#ifndef HERWIG_TauTo2MesonAnalysis_H
#define HERWIG_TauTo2MesonAnalysis_H
//
// This is the declaration of the TauTo2MesonAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The TauTo2MesonAnalysis class is designed to perform the analysis of the 
 * mass distribution of the hadronic decay products of the \f$\tau\f$ in the decays
 * \f$\tau^\pm\to\nu_\tau\{\pi^\pm\pi^0,K^\pm\pi^0,K^0\pi^\pm,K^\pm\eta,K^\pm K^0\}\f$.
 * In order to work the \f$\pi^0\f$, \f$K^0\f$, \f$K^\pm\f$ and \f$\eta\f$ should be
 * set stable.
 *
 * The mass spectrum of the \f$pi^\pm\pi^0\f$ final state is compared to CLEO and
 * BELLE data.
 *
 * @see \ref TauTo2MesonAnalysisInterfaces "The interfaces"
 * defined for TauTo2MesonAnalysis.
 */
class TauTo2MesonAnalysis: public AnalysisHandler {

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
   * Indicates that this is a concrete class.
   */
  static NoPIOClassDescription<TauTo2MesonAnalysis> initTauTo2MesonAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TauTo2MesonAnalysis & operator=(const TauTo2MesonAnalysis &);

private:

  /**
   *  Histograms for the mass plots
   */
  //@{
  /**
   *  Mass of the pions in \f$\tau\pm\to\nu_\tau\pi^0\pi^\pm\f$ compared to BELLE
   * and CLEO data.
   */
  HistogramPtr _m2pipiBELLE,_mpipiCLEO;

  /**
   *  Mass of the Kaons and pions in \f$\tau\to K\pi\f$
   */
  HistogramPtr _m2KpiA,_m2KpiB,_mKpiA,_mKpiB,_m2KpiC,_m2KpiD;

  /**
   *  Mass of the \f$K\eta\f$
   */
  HistogramPtr _m2Keta,_mKeta;

  /**
   *  Mass of the \f$KK\f$
   */
  HistogramPtr _m2KK,_mKK;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TauTo2MesonAnalysis. */
template <>
struct BaseClassTrait<Herwig::TauTo2MesonAnalysis,1> {
  /** Typedef of the first base class of TauTo2MesonAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TauTo2MesonAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TauTo2MesonAnalysis>
  : public ClassTraitsBase<Herwig::TauTo2MesonAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TauTo2MesonAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TauTo2MesonAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class TauTo2MesonAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwTauAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_TauTo2MesonAnalysis_H */
