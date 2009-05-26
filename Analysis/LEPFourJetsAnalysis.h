// -*- C++ -*-
//
// LEPFourJetsAnalysis.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
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
#include "Herwig++/Interfaces/KtJetInterface.h"
#include "KtJet/KtJet.h"
#include "KtJet/KtLorentzVector.h"
#include "Herwig++/Utilities/Histogram.h"
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
  inline double cosChiBZ(vector<Lorentz5Momentum>);

  /**
   *  Compute \f$\cos\Phi_{KSW}\f$.
   */ 
  inline double cosPhiKSW(vector<Lorentz5Momentum>);
  
  /**
   *  Compute \f$\cos\Theta_{NR}\f$
   */
  inline double cosThetaNR(vector<Lorentz5Momentum>); 

  /**
   *  Compute \f$\cos\alpha_{34}\f$
   */
  inline double cosAlpha34(vector<Lorentz5Momentum>); 
  //@}

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
  static NoPIOClassDescription<LEPFourJetsAnalysis> initLEPFourJetsAnalysis;

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
   *  The interface between Herwig++ and KtJet
   */
  Herwig::KtJetInterface _kint;
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
  static string library() { return "HwKtJet.so HwLEPJetAnalysis.so"; }
};

/** @endcond */

}

#include "LEPFourJetsAnalysis.icc"

#endif /* HERWIG_LEPFourJetsAnalysis_H */
