// -*- C++ -*-
//
// SingleParticleAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SingleParticleAnalysis_H
#define HERWIG_SingleParticleAnalysis_H
//
// This is the declaration of the SingleParticleAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Utilities/Histogram.h"
#include "EventShapes.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The SingleParticleAnalysis class performs the analysis for
 * single particle variables and is intended as a slave handler for the 
 *  EventShapesMasterAnalysis class.
 *
 * @see \ref SingleParticleAnalysisInterfaces "The interfaces"
 * defined for SingleParticleAnalysis.
 */
class SingleParticleAnalysis: public AnalysisHandler {

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);
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
  static ClassDescription<SingleParticleAnalysis> initSingleParticleAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SingleParticleAnalysis & operator=(const SingleParticleAnalysis &);

private:

  /**
   * Histogram for the rapidity distribution with respect to the thrust axis
   */
  HistogramPtr _yT;

  /**
   * Histogram for the rapidity distribution with respect to the sphericity axis
   */
  HistogramPtr _yS;

  /**
   * Histogram for the \f$p_{T,in}\f$ distribution with respect to the thrust axis
   */
  HistogramPtr _ptinT;

  /**
   * Histogram for the \f$p_{T,out}\f$ distribution with respect to the thrust axis
   */
  HistogramPtr _ptoutT;

  /**
   * Histogram for the \f$p_{S,in}\f$ distribution with respect to the sphericity axis
   */
  HistogramPtr _ptinS;

  /**
   * Histogram for the \f$p_{S,out}\f$ distribution with respect to the sphericity axis
   */
  HistogramPtr _ptoutS;

  /**
   * Histogram for the number of charged particles
   */
  HistogramPtr _nch;

  /**
   *  Object which calculates the event shapes
   */
  EventShapesPtr _shapes;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SingleParticleAnalysis. */
template <>
struct BaseClassTrait<Herwig::SingleParticleAnalysis,1> {
  /** Typedef of the first base class of SingleParticleAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SingleParticleAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SingleParticleAnalysis>
  : public ClassTraitsBase<Herwig::SingleParticleAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SingleParticleAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SingleParticleAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so HwLEPAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SingleParticleAnalysis_H */
