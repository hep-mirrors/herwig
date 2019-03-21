// -*- C++ -*-
#ifndef HERWIG_DDalitzAnalysis_H
#define HERWIG_DDalitzAnalysis_H
//
// This is the declaration of the DDalitzAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/EventRecord/Particle.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the DDalitzAnalysis class.
 *
 * @see \ref DDalitzAnalysisInterfaces "The interfaces"
 * defined for DDalitzAnalysis.
 */
class DDalitzAnalysis: public AnalysisHandler {

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
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

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

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Find the stable decay products
   */
  void findChildren(tPPtr,ParticleVector &);

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DDalitzAnalysis> initDDalitzAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DDalitzAnalysis & operator=(const DDalitzAnalysis &) = delete;

private:

  /**
   *  Histograms for \f$D^0\to \bar{K}^0\pi^+\pi^-\f$
   */
  //@{
  /**
   * \f$m^2_+\f$
   */
  HistogramPtr _m2plus1;

  /**
   * \f$m^2_+\f$
   */
  HistogramPtr _m2minus1;

  /**
   * \f$m^2_{\pi\pi}\f$
   */
  HistogramPtr _m2pipi1;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points1;
  //@}

  /**
   *  Histograms for \f$D^0\to K^-\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Histogram for the \f$K^-\pi^+\f$ mass
   */
  HistogramPtr _m2minus2;

  /**
   *  Histogram for the \f$\pi^+\pi^0\f$ mass
   */
  HistogramPtr _m2pipi2;

  /**
   *  Histogram for the \f$K^-\pi^0\f$ mass
   */
  HistogramPtr _m2neutral2;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points2;
  //@}

  /**
   *  Histograms for \f$D^+\to K^-\pi^+\pi^+\f$
   */
  //@{
  /**
   *  Histogram for \f$K^-\pi^+\f$ low
   */
  HistogramPtr _mKpilow3;

  /**
   *  Histogram for \f$K^-\pi^+\f$ high
   */
  HistogramPtr _mKpihigh3;

  /**
   *  Histogram for \f$K^-\pi^+\f$ all
   */
  HistogramPtr _mKpiall3;

  /**
   *  Histogram for \f$\pi^+\pi^-\f$
   */
  HistogramPtr _mpipi3;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points3;
  //@}

  /**
   *  Histograms for \f$D^+\to\bar{K}^0\pi^+\pi^0\f$
   */
  //@{
  /**
   *  Histogram for the \f$\bar{K}^0\pi^+\f$ mass
   */
  HistogramPtr _m2Kpip4;

  /**
   *  Histogram for the \f$\pi^+\pi^0\f$ mass
   */
  HistogramPtr _m2pipi4;

  /**
   *  Histogram for the \f$\bar{K}^0\pi^0\f$ mass
   */
  HistogramPtr _m2Kpi04;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points4;
  //@}

  /**
   *  Histograms for \f$D^+\to K^+\pi^-\pi^+\f$
   */
  //@{
  /**
   *  Histogram for \f$K^+\pi^-\f$
   */
  HistogramPtr _mkppim5;

  /**
   *  Histogram for \f$K^+\pi^+\f$
   */
  HistogramPtr _mkppip5;

  /**
   *  Histogram for \f$\pi^+\pi^-\f$
   */
  HistogramPtr _mpippim5;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points5;
  //@}

  /**
   *  Histograms for \f$D_s^+\to K^+\pi^-\pi^+\f$
   */
  //@{
  /**
   *  Histogram for \f$K^+\pi^-\f$
   */
  HistogramPtr _mkppim6;

  /**
   *  Histogram for \f$K^+\pi^+\f$
   */
  HistogramPtr _mkppip6;

  /**
   *  Histogram for \f$\pi^+\pi^-\f$
   */
  HistogramPtr _mpippim6;

  /**
   *  Vectors for the Dalitz plot
   */
  vector<pair<Energy2,Energy2> > _points6;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DDalitzAnalysis. */
template <>
struct BaseClassTrait<Herwig::DDalitzAnalysis,1> {
  /** Typedef of the first base class of DDalitzAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DDalitzAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DDalitzAnalysis>
  : public ClassTraitsBase<Herwig::DDalitzAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DDalitzAnalysis"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DDalitzAnalysis is implemented. It may also include several, space-separated,
   * libraries if the class DDalitzAnalysis depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDecayAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DDalitzAnalysis_H */
