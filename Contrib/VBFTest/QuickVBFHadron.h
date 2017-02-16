// -*- C++ -*-
#ifndef HERWIG_QuickVBFHadron_H
#define HERWIG_QuickVBFHadron_H
//
// This is the declaration of the QuickVBFHadron class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the QuickVBFHadron class.
 *
 * @see \ref QuickVBFHadronInterfaces "The interfaces"
 * defined for QuickVBFHadron.
 */
class QuickVBFHadron: public AnalysisHandler {

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

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
  static NoPIOClassDescription<QuickVBFHadron> initQuickVBFHadron;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QuickVBFHadron & operator=(const QuickVBFHadron &);

private:

  HistogramPtr _mH  ,_yH  ,_pTH[2]  ,_phiH  ;
  HistogramPtr _yjet,_pTjet[2],_phijet;
  HistogramPtr _mjj;

  static bool _doOnce;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QuickVBFHadron. */
template <>
struct BaseClassTrait<Herwig::QuickVBFHadron,1> {
  /** Typedef of the first base class of QuickVBFHadron. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QuickVBFHadron class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QuickVBFHadron>
  : public ClassTraitsBase<Herwig::QuickVBFHadron> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QuickVBFHadron"; }
  /**
   * The name of a file containing the dynamic library where the class
   * QuickVBFHadron is implemented. It may also include several, space-separated,
   * libraries if the class QuickVBFHadron depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "VBFAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_QuickVBFHadron_H */
