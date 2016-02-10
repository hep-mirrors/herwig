// -*- C++ -*-
//
// LEPFourJetsAnalysis.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_LEPFourJetsAnalysis_H
#define HERWIG_LEPFourJetsAnalysis_H
//
// This is the declaration of the LEPFourJetsAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "Herwig/Utilities/Histogram.h"
#include "ThePEG/Repository/CurrentGenerator.h"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Analysis
 * The LEPFourJetsAnalysis class performs analysis for four jet angles and
 * compares them to LEP data.
 *
 * @see \ref LEPFourJetsAnalysisInterfaces "The interfaces"
 * defined for LEPFourJetsAnalysis.
 */
class LEPFourJetsAnalysis: public AnalysisHandler {

public:

  /**
   * Default constructor
   */
  LEPFourJetsAnalysis ()
    : _ca34(), _cchiBZ(), _cphiKSW(), _cthNR(),
      _charged(true) {}

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
   *  Methods to compute the four jet angles, assumes the jets are energy ordered
   */
  //@{
  /**
   *  Compute \f$\cos\chi_{BZ}\f$
   */
  double cosChiBZ(vector<Lorentz5Momentum> p) {
    if (p.size() == 4) {
      ThreeVector<Energy2> v1 = p[0].vect().cross(p[1].vect());
      ThreeVector<Energy2> v2 = p[2].vect().cross(p[3].vect());
      return cos(v1.angle(v2)); 
    } 
    else return 123;
  }

  /**
   *  Compute \f$\cos\Phi_{KSW}\f$.
   */ 
  double cosPhiKSW(vector<Lorentz5Momentum> p) {
    if (p.size() == 4) {
      ThreeVector<Energy2> v1 = p[0].vect().cross(p[3].vect());
      ThreeVector<Energy2> v2 = p[1].vect().cross(p[2].vect());
      double alpha1 = v1.angle(v2);
      v1 = p[0].vect().cross(p[2].vect());
      v2 = p[1].vect().cross(p[3].vect());
      double alpha2 = v1.angle(v2);
      return cos((alpha1+alpha2)/2.);
    } 
    else return 123; 
  }

  /**
   *  Compute \f$\cos\Theta_{NR}\f$
   */
  double cosThetaNR(vector<Lorentz5Momentum> p) {
    if (p.size() == 4) {
      ThreeVector<Energy> v1 = p[0].vect() - p[1].vect();
      ThreeVector<Energy> v2 = p[2].vect() - p[3].vect();
      return cos(v1.angle(v2));
    }
    else return 123; 
  }

  /**
   *  Compute \f$\cos\alpha_{34}\f$
   */
  double cosAlpha34(std::vector<Lorentz5Momentum> p) {
    if (p.size() == 4)
      return cos(p[2].vect().angle(p[3].vect()));
    else 
      return 123; 
  }
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
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
  static ClassDescription<LEPFourJetsAnalysis> initLEPFourJetsAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  LEPFourJetsAnalysis & operator=(const LEPFourJetsAnalysis &);

private:

  /**
   * Histogram for the \f$\cos\alpha_{34}\f$ distribution
   */
  HistogramPtr _ca34;

  /**
   * Histogram for the \f$\cos\chi_{BZ}\f$ distribution
   */
  HistogramPtr _cchiBZ;

  /**
   * Histogram for the \f$\cos\Phi_{KSW}\f$ distribution
   */
  HistogramPtr _cphiKSW;

  /**
   * Histogram for the \f$\cos\Theta_{NR}\f$ distribution
   */
  HistogramPtr _cthNR;

  /**
   * Use charged particles only
   */
  bool _charged;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of LEPFourJetsAnalysis. */
template <>
struct BaseClassTrait<Herwig::LEPFourJetsAnalysis,1> {
  /** Typedef of the first base class of LEPFourJetsAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the LEPFourJetsAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::LEPFourJetsAnalysis>
  : public ClassTraitsBase<Herwig::LEPFourJetsAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::LEPFourJetsAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the LEPFourJetsAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwLEPJetAnalysis.so"; }
};

/** @endcond */

}

#endif /* HERWIG_LEPFourJetsAnalysis_H */
