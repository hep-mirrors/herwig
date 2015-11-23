// -*- C++ -*-
#ifndef HERWIG_EtaDecayAnalysis_H
#define HERWIG_EtaDecayAnalysis_H
//
// This is the declaration of the EtaDecayAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The EtaDecayAnalysis class performs the analysis of the decays of 
 * the \f$\eta\f$ and \f$\eta'\f$ mesons.
 *
 * @see \ref EtaDecayAnalysisInterfaces "The interfaces"
 * defined for EtaDecayAnalysis.
 */
class EtaDecayAnalysis: public AnalysisHandler {

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
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<EtaDecayAnalysis> initEtaDecayAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EtaDecayAnalysis & operator=(const EtaDecayAnalysis &);

private:

  /**
   *  Histograms for the decay \f$\eta\to\pi^0\gamma\gamma\f$
   */
  //@{
  /**
   * Histogram for the mass of \f$\gamma\gamma\f$
   */
  vector<HistogramPtr> _mgammagamma;

  /**
   * Histogrma for the mass of \f$\pi^0\gamma\f$
   */
  vector<HistogramPtr> _mpi0gamma;
  //@}

  /**
   *  Histograms for the decay \f$\eta\to\pi^+\pi^-\gamma\f$
   */
  //@{
  /**
   *  Histogram for the mass of \f$\pi^+\gamma\f$
   */
  vector<HistogramPtr> _mpipgamma;

  /**
   *  Histogram for the mass of \f$\pi^-\gamma\f$
   */
  vector<HistogramPtr> _mpimgamma;

  /**
   *  Histogram for the mass of \f$\pi^+\pi^-\f$
   */
  vector<HistogramPtr> _mpippim;

  /**
   *  Histogram for the photon energy
   */
  vector<HistogramPtr> _photonenergy;
  //@}

  /**
   * Histograms for the decay \f$\eta\pi\pi\pi\f$ and \f$\eta'\to\eta\pi\pi\f$.
   */
  //@{
  /**
   *  Histogram for the mass of \f$\pi^+\pi^-\f$
   */
  vector<HistogramPtr> _dpippim;

  /**
   *  Histogram for the mass of \f$\pi^0\pi^0\f$
   */
  vector<HistogramPtr> _dpi0pi0;

  /**
   *  Histogram for the mass of \f$\pi^0\pi^+\f$
   */
  vector<HistogramPtr> _dpi0pip;

  /**
   *  Histogram for the mass of \f$\pi^0\pi^-\f$
   */
  vector<HistogramPtr> _dpi0pim;

  /**
   *  Histogram for the mass of \f$\pi^+\eta\f$
   */
  vector<HistogramPtr> _dpipeta;

  /**
   *  Histogram for the mass of \f$\pi^-\eta\f$
   */
  vector<HistogramPtr> _dpimeta;

  /**
   *  Histogram for the mass of \f$\pi^0\eta\f$
   */
  vector<HistogramPtr> _dpi0eta;

  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EtaDecayAnalysis. */
template <>
struct BaseClassTrait<Herwig::EtaDecayAnalysis,1> {
  /** Typedef of the first base class of EtaDecayAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EtaDecayAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EtaDecayAnalysis>
  : public ClassTraitsBase<Herwig::EtaDecayAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EtaDecayAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * EtaDecayAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class EtaDecayAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_EtaDecayAnalysis_H */
