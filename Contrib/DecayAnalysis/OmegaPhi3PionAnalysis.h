// -*- C++ -*-
#ifndef HERWIG_OmegaPhi3PionAnalysis_H
#define HERWIG_OmegaPhi3PionAnalysis_H
//
// This is the declaration of the OmegaPhi3PionAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the OmegaPhi3PionAnalysis class.
 *
 * @see \ref OmegaPhi3PionAnalysisInterfaces "The interfaces"
 * defined for OmegaPhi3PionAnalysis.
 */
class OmegaPhi3PionAnalysis: public AnalysisHandler {

public:

  /**
   *  Default Constructor
   */
  inline OmegaPhi3PionAnalysis() : _nmax(50000) {}

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<OmegaPhi3PionAnalysis> initOmegaPhi3PionAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  OmegaPhi3PionAnalysis & operator=(const OmegaPhi3PionAnalysis &) = delete;

private:

  /**
   *  Histogram for the \f$x\f$-values
   */
  vector<HistogramPtr> _xhist;

  /**
   *  Histogram for the \f$y\f$-values
   */
  vector<HistogramPtr> _yhist;

  /**
   *  Histograms for the masses
   */
  //@{
  /**
   *  The mass of the \f$\rho^+\f$
   */
  vector<HistogramPtr> _mplus;

  /**
   *  The mass of the \f$\rho^-\f$
   */
  vector<HistogramPtr> _mminus;

  /**
   *  The mass of the \f$\rho^0\f$
   */
  vector<HistogramPtr> _m0;
  //@}
  /**
   *  Vectors to store the \f$x\f$ and\f$y\f$ values for
   */
  //@{
  /**
   *  The \f$x\f$ value
   */
  vector<vector<Energy> > _xvalue;

  /**
   *  The \f$y\f$ value
   */
  vector<vector<Energy> > _yvalue;
  //@}

  /**
   *  Maximum number of points for the Dalitz plots
   */
  unsigned int _nmax;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of OmegaPhi3PionAnalysis. */
template <>
struct BaseClassTrait<Herwig::OmegaPhi3PionAnalysis,1> {
  /** Typedef of the first base class of OmegaPhi3PionAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the OmegaPhi3PionAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::OmegaPhi3PionAnalysis>
  : public ClassTraitsBase<Herwig::OmegaPhi3PionAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::OmegaPhi3PionAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * OmegaPhi3PionAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class OmegaPhi3PionAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_OmegaPhi3PionAnalysis_H */
