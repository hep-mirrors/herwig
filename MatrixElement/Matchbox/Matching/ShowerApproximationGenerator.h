// -*- C++ -*-
//
// ShowerApproximationGenerator.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_ShowerApproximationGenerator_H
#define Herwig_ShowerApproximationGenerator_H
//
// This is the declaration of the ShowerApproximationGenerator class.
//

#include "ThePEG/Handlers/StepHandler.h"
#include "Herwig/MatrixElement/Matchbox/Matching/ShowerApproximationKernel.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/MatchboxPhasespace.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief ShowerApproximationGenerator generates emissions according to a
 * shower approximation entering a NLO matching.
 *
 */
class ShowerApproximationGenerator: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ShowerApproximationGenerator();

  /**
   * The destructor.
   */
  virtual ~ShowerApproximationGenerator();
  //@}

public:

  /** @name Virtual functions required by the StepHandler class. */
  //@{
  /**
    * The main function called by the EventHandler class to
    * perform a step. Given the current state of an Event, this function
    * performs the event generation step and includes the result in a new
    * Step object int the Event record.
    * @param eh the EventHandler in charge of the Event generation.
    * @param tagged if not empty these are the only particles which should
    * be considered by the StepHandler.
    * @param hint a Hint object with possible information from previously
    * performed steps.
    * @throws Veto if the StepHandler requires the current step to be discarded.
    * @throws Stop if the generation of the current Event should be stopped
    * after this call.
    * @throws Exception if something goes wrong.
    */
  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint);
  //@}

public:

  /**
   * Fill information on the Born process
   */
  bool prepare(bool didproject);

  /**
   * Generate a Born phase space point while kernels are being
   * presampled
   */
  bool generate(const vector<double>&);

  /**
   * Restore information on the Born process
   */
  void restore();

protected:

  /**
   * Generate a momentum fraction for the given parton species
   */
  double generateFraction(tcPDPtr, double, double) const;

  /**
   * Invert the momentum fraction for the given parton species
   */
  double invertFraction(tcPDPtr, double, double) const;

  /**
   * Return the pt cut to be applied for final-final dipoles.
   */
  Energy ffPtCut() const { return theShowerApproximation->ffPtCut(); }

  /**
   * Return the pt cut to be applied for final-initial dipoles.
   */
  Energy fiPtCut() const { return theShowerApproximation->fiPtCut(); }

  /**
   * Return the pt cut to be applied for initial-initial dipoles.
   */
  Energy iiPtCut() const { return theShowerApproximation->iiPtCut(); }

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The shower approximation to consider
   */
  Ptr<ShowerApproximation>::ptr theShowerApproximation;

  /**
   * The (invertible) phase space generator to use
   */
  Ptr<MatchboxPhasespace>::ptr thePhasespace;

  /**
   * Map hard processes to the respective kernels.
   */
  map<cPDVector,set<Ptr<ShowerApproximationKernel>::ptr> > theKernelMap;

  /**
   * The number of points to presample this splitting generator.
   */
  unsigned long thePresamplingPoints;

  /**
   * The maximum number of trials to generate a splitting.
   */
  unsigned long theMaxTry;

  /**
   * Return the number of accepted points after which the grid should
   * be frozen
   */
  unsigned long theFreezeGrid;

  /**
   * The last external Born XComb dealt with
   */
  tStdXCombPtr lastIncomingXComb;

  // the next three are filled from the incoming xcomb in the prepare method

  /**
   * The last internal Born matrix element dealt with
   */
  Ptr<MatchboxMEBase>::ptr theLastBornME;

  /**
   * The last Born phase space point
   */
  vector<Lorentz5Momentum> theLastMomenta;

  /**
   * The last Born phase space point used while presampling
   */
  vector<Lorentz5Momentum> theLastPresamplingMomenta;

  /**
   * The random numbers which have produced the last Born phase space
   * point.
   */
  vector<double> theLastRandomNumbers;

  /**
   * The last internal Born XComb dealt with
   */
  tStdXCombPtr theLastBornXComb;

  /**
   * The last internal incoming partons dealt with
   */
  PPair theLastPartons;

  /**
   * True, if sampler should apply compensation
   */
  bool theDoCompensate;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ShowerApproximationGenerator & operator=(const ShowerApproximationGenerator &) = delete;

};

}

#endif /* Herwig_ShowerApproximationGenerator_H */
