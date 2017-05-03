// -*- C++ -*-
//
// BasicConsistency.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_BasicConsistency_H
#define THEPEG_BasicConsistency_H
//
// This is the declaration of the BasicConsistency class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Analysis
 * The BasicConsistency class is a simple analysis which performs a basic
 * analysis of the event checking that energy, momentum and charge are
 * conserved and no quarks or clusters are final-state particles.
 *
 * @see \ref BasicConsistencyInterfaces "The interfaces"
 * defined for BasicConsistency.
 */
class BasicConsistency: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  BasicConsistency();

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BasicConsistency> initBasicConsistency;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BasicConsistency & operator=(const BasicConsistency &);

private:

  /**
   *  Maximum momentum deviation
   */
  Energy _epsmom;

  /**
   *  check for quarks
   */
  bool _checkquark;

  /**
   *  check for charge conservation
   */
  bool _checkcharge;

  /**
   *  Check for clusters in  the final-state
   */
  bool _checkcluster;

  /**
   *  Check the branching ratios
   */
  bool _checkBR;

  /**
   *  Maximum absolute momentum deviation before warning
   */
  Energy _absolutemomentumtolerance;

  /**
   *  Maximum momentum deviation relative to beam energy before warning
   */
  double _relativemomentumtolerance;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BasicConsistency. */
template <>
struct BaseClassTrait<Herwig::BasicConsistency,1> {
  /** Typedef of the first base class of BasicConsistency. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BasicConsistency class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BasicConsistency>
  : public ClassTraitsBase<Herwig::BasicConsistency> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BasicConsistency"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the BasicConsistency class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_BasicConsistency_H */
