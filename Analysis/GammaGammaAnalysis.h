// -*- C++ -*-
#ifndef HERWIG_GammaGammaAnalysis_H
#define HERWIG_GammaGammaAnalysis_H
//
// This is the declaration of the GammaGammaAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "GammaGammaAnalysis.fh"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * GammaGammaAnalysis is for the analysis of events with a pair of hard
 * photons produced.  These are selected as the two highest pt photons
 * in the final state of the event.  A topdrawer file with histograms
 * is written to the working directory.  
 *
 * @see \ref GammaGammaAnalysisInterfaces "The interfaces"
 * defined for GammaGammaAnalysis.
 */
class GammaGammaAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline GammaGammaAnalysis();

  /**
   * The destructor.
   */
  virtual ~GammaGammaAnalysis();
  //@}

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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GammaGammaAnalysis> initGammaGammaAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GammaGammaAnalysis & operator=(const GammaGammaAnalysis &);

private:
  /**
   *   \f$p_T\f$ of the photon
   */
  Histogram _ptharder;
  Histogram _ptsofter;
  Histogram _ptpair;

  /**
   *   Energy of the photons
   */
  Histogram _Eharder;
  Histogram _Esofter;
  Histogram _Epair;

  /**
   *  Rapidity of the photons
   */
  Histogram _rapharder;
  Histogram _rapsofter;
  Histogram _rappair;

  /**
   *  Azimuth of the photons
   */
  Histogram _phiharder;
  Histogram _phisofter;
  Histogram _deltaphi;
  
  /**
   *  Invariant mass of the pair
   */
  Histogram _mpair;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GammaGammaAnalysis. */
template <>
struct BaseClassTrait<Herwig::GammaGammaAnalysis,1> {
  /** Typedef of the first base class of GammaGammaAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GammaGammaAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GammaGammaAnalysis>
  : public ClassTraitsBase<Herwig::GammaGammaAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::GammaGammaAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GammaGammaAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#include "GammaGammaAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GammaGammaAnalysis.tcc"
#endif

#endif /* HERWIG_GammaGammaAnalysis_H */
