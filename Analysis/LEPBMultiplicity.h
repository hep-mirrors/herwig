// -*- C++ -*-
//
// LEPBMultiplicity.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LEPBMultiplicity_H
#define HERWIG_LEPBMultiplicity_H
//
// This is the declaration of the LEPBMultiplicity class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace Herwig {

using namespace ThePEG;


/** \ingroup Analysis
 *  Struct for the multiplcity data
 */
struct BranchingInfo {
  /**
   *  Default constructor
   * @param mult  The observed multiplcity.
   * @param error The error on the observed multiplicity
   */
  BranchingInfo(double mult=0.,double error=0.);

  /**
   *  The observed multiplicity
   */
  double obsBranching;

  /**
   *  The error on the observed multiplicity
   */
  double obsError;

  /**
   *  Number of particles of this type
   */
  long actualCount;

  /**
   *  Sum of squares of number per event for error
   */
  double sumofsquares;

  /**
   *  The average fraction per quark
   * @param N The number of events
   * @param den The denominator to give the fraction
   */
  double simBranching(long N,BranchingInfo den=BranchingInfo());

  /**
   *  The error on the average number per event
   * @param N The number of events 
   * @param den The denominator to give the fraction
   */
  double simError(long N,BranchingInfo den=BranchingInfo());

  /**
   * Is the result more than \f$3\sigma\f$ from the experimental result
   * @param N The number of events
   * @param den The denominator to give the fraction
   */
  double nSigma(long N,BranchingInfo den=BranchingInfo());

  /**
   * Plot standard error in a simple barchart
   * @param N The number of events
   * @param den The denominator to give the fraction
   */
  string bargraph(long N,BranchingInfo den=BranchingInfo());
};

/**
 * The LEPBBMultiplicity class is designed to compare the production
 * rates of \f$B^+\f$, \f$B^0\f$, \f$B^0_s\f$ and B-baryons in
 * B events at LEP
 *
 * @see \ref LEPBMultiplicityInterfaces "The interfaces"
 * defined for LEPBMultiplicity.
 */
class LEPBMultiplicity: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  LEPBMultiplicity();

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
  static NoPIOClassDescription<LEPBMultiplicity> initLEPBMultiplicity;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LEPBMultiplicity & operator=(const LEPBMultiplicity &);

private:

  /**
   *  Map of PDG codes to multiplicity info
   */
  map<long,BranchingInfo> _data;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LEPBMultiplicity. */
template <>
struct BaseClassTrait<Herwig::LEPBMultiplicity,1> {
  /** Typedef of the first base class of LEPBMultiplicity. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LEPBMultiplicity class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LEPBMultiplicity>
  : public ClassTraitsBase<Herwig::LEPBMultiplicity> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LEPBMultiplicity"; }
  /**
   * The name of a file containing the dynamic library where the class
   * LEPBMultiplicity is implemented. It may also include several, space-separated,
   * libraries if the class LEPBMultiplicity depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LEPBMultiplicity_H */
