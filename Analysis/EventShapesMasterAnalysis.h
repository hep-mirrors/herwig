// -*- C++ -*-
//
// EventShapesMasterAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_EventShapesMasterAnalysis_H
#define HERWIG_EventShapesMasterAnalysis_H
//
// This is the declaration of the EventShapesMasterAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "EventShapes.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The EventShapesMasterAnalysis class is responsible for extracting the 
 * final state particles from the event record and setting up the object which
 * calculates the event shapes. This is done so that the EventShapes class which 
 * contains a lot of slow code isn't rerun unnecessarily.
 *
 * @see \ref EventShapesMasterAnalysisInterfaces "The interfaces"
 * defined for EventShapesMasterAnalysis.
 */
class EventShapesMasterAnalysis: public AnalysisHandler {

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<EventShapesMasterAnalysis> initEventShapesMasterAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EventShapesMasterAnalysis & operator=(const EventShapesMasterAnalysis &);

private:

  /**
   *  Pointer to the EventShapes object
   */
  EventShapesPtr _shapes;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of EventShapesMasterAnalysis. */
template <>
struct BaseClassTrait<Herwig::EventShapesMasterAnalysis,1> {
  /** Typedef of the first base class of EventShapesMasterAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the EventShapesMasterAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::EventShapesMasterAnalysis>
  : public ClassTraitsBase<Herwig::EventShapesMasterAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::EventShapesMasterAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the EventShapesMasterAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_EventShapesMasterAnalysis_H */
