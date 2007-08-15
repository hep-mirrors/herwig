// -*- C++ -*-
#ifndef THEPEG_BasicConsistency_H
#define THEPEG_BasicConsistency_H
//
// This is the declaration of the BasicConsistency class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "BasicConsistency.fh"

namespace Herwig {
  using namespace ThePEG;

/**
 * Here is the documentation of the BasicConsistency class.
 *
 * @see \ref BasicConsistencyInterfaces "The interfaces"
 * defined for BasicConsistency.
 */
class BasicConsistency: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  inline BasicConsistency();

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

#include "BasicConsistency.icc"

#endif /* THEPEG_BasicConsistency_H */
