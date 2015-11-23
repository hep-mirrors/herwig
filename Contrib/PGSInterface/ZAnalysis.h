// -*- C++ -*-
#ifndef HERWIG_ZAnalysis_H
#define HERWIG_ZAnalysis_H
//
// This is the declaration of the ZAnalysis class.
//

#include "PGSInterface.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The ZAnalysis class is a simple example of using the interface
 * to PGS to look at the Z mass at the hadron and detector levels
 *
 * @see \ref ZAnalysisInterfaces "The interfaces"
 * defined for ZAnalysis.
 */
class ZAnalysis: public PGSInterface {

public:

  /**
   * The default constructor.
   */
  ZAnalysis();

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
  static NoPIOClassDescription<ZAnalysis> initZAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ZAnalysis & operator=(const ZAnalysis &);

private:

  /**
   *  Z mass at the hadron level
   */
  Histogram ZmassHadron_;

  /**
   *  Z mass at the detector level
   */
  Histogram ZmassDetector_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ZAnalysis. */
template <>
struct BaseClassTrait<Herwig::ZAnalysis,1> {
  /** Typedef of the first base class of ZAnalysis. */
  typedef Herwig::PGSInterface NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ZAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::ZAnalysis>
  : public ClassTraitsBase<Herwig::ZAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::ZAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ZAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class ZAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwPGSInterface.so ZAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_ZAnalysis_H */
