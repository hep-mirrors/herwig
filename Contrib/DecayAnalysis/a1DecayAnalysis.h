// -*- C++ -*-
#ifndef HERWIG_a1DecayAnalysis_H
#define HERWIG_a1DecayAnalysis_H
//
// This is the declaration of the a1DecayAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the a1DecayAnalysis class.
 *
 * @see \ref a1DecayAnalysisInterfaces "The interfaces"
 * defined for a1DecayAnalysis.
 */
class a1DecayAnalysis: public AnalysisHandler {

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
   *  find the pions produced in the decay
   */
  void findPions(tPPtr part,ParticleVector & pions);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<a1DecayAnalysis> inita1DecayAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  a1DecayAnalysis & operator=(const a1DecayAnalysis &) = delete;

private:

  /**
   *  Histograms for \f$a_1^0\to\pi^0\pi^0\pi^0\f$
   */
  HistogramPtr _hist0;

  /**
   *  Histograms for \f$a_1^+\to\pi^0\pi^0\pi^+\f$
   */
  //@{
  /**
   * Mass of the \f$\pi^0\pi^0\f$ pair 
   */
  HistogramPtr _hist1A;

  /**
   * Mass of the \f$\pi^0\pi^+\f$ pair
   */
  HistogramPtr _hist1B;
  //@}

  /**
   *  Histograms for \f$a_1^0\to\pi^+\pi^-\pi^0\f$
   */
  //@{
  /**
   * Mass of the \f$\pi^+\pi^-\f$ pair 
   */
  HistogramPtr _hist2A;

  /**
   * Mass of the \f$\pi^+\pi^0\f$ pair
   */
  HistogramPtr _hist2B;

  /**
   * Mass of the \f$\pi^-\pi^0\f$ pair
   */
  HistogramPtr _hist2C;
  //@}

  /**
   *  Histograms for \f$a_1^+\to\pi^+\pi^+\pi^-\f$
   */
  //@{
  /**
   * Mass of the \f$\pi^+\pi^+\f$ pair 
   */
  HistogramPtr _hist3A;

  /**
   * Mass of the \f$\pi^+\pi^-\f$ pair
   */
  HistogramPtr _hist3B;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of a1DecayAnalysis. */
template <>
struct BaseClassTrait<Herwig::a1DecayAnalysis,1> {
  /** Typedef of the first base class of a1DecayAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the a1DecayAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::a1DecayAnalysis>
  : public ClassTraitsBase<Herwig::a1DecayAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::a1DecayAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * a1DecayAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class a1DecayAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_a1DecayAnalysis_H */
