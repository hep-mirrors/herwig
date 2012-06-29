// -*- C++ -*-
//
// PowhegSplittingGenerator.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PowhegSplittingGenerator_H
#define HERWIG_PowhegSplittingGenerator_H
//
// This is the declaration of the PowhegSplittingGenerator class.
//

#include "ThePEG/Handlers/StepHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/MatrixElement/Matchbox/Powheg/PowhegSplittingKernel.h"
#include "Herwig++/Exsample2/exsample/exponential_generator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief PowhegSplittingGenerator generates the POWHEG
 * real emission contribution.
 *
 * @see \ref PowhegSplittingGeneratorInterfaces "The interfaces"
 * defined for PowhegSplittingGenerator.
 */
class PowhegSplittingGenerator: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  PowhegSplittingGenerator();

  /**
   * The destructor.
   */
  virtual ~PowhegSplittingGenerator();
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
   * Generate a splitting given the current event handler
   */
  bool generate(EventHandler & eh);

  /**
   * Return the last selected splitting
   */
  Ptr<PowhegSplittingKernel>::tcptr lastSplitting() const {
    return theLastSplitting;
  }

  /**
   * Return the last selected splitting
   */
  Ptr<PowhegSplittingKernel>::tptr lastSplitting() {
    return theLastSplitting;
  }

  /**
   * Indicate that the next event should
   * be discarded owing to presampling.
   */
  void setDiscardNext(bool f = true) { theDiscardNext = f; }

  /**
   * Return true, if the next event should
   * be discarded owing to presampling.
   */
  bool discardNext() const { return theDiscardNext; }

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
   * The pt cut to be applied for final-final dipoles.
   */
  Energy theFFPtCut;

  /**
   * An optional screening scale for final-final dipoles; see
   * DipoleSplittingKernel
   */
  Energy theFFScreeningScale;

  /**
   * The pt cut to be applied for final-initial dipoles.
   */
  Energy theFIPtCut;

  /**
   * An optional screening scale for final-initial dipoles; see
   * DipoleSplittingKernel
   */
  Energy theFIScreeningScale;

  /**
   * The pt cut to be applied for initial-initial dipoles.
   */
  Energy theIIPtCut;

  /**
   * An optional screening scale for initial-initial dipoles; see
   * DipoleSplittingKernel
   */
  Energy theIIScreeningScale;

  /**
   * Discard events which did not radiate.
   */
  bool discardNoEmissions;

  typedef
  exsample::exponential_generator<PowhegSplittingKernel,UseRandom>
  ExponentialGenerator;

  typedef
  exsample::exponential_generator<PowhegSplittingKernel,UseRandom>*
  ExponentialGeneratorPtr;

  typedef
  multimap<XCombPtr,pair<Ptr<PowhegSplittingKernel>::ptr,ExponentialGeneratorPtr> >
  GeneratorMap;

  /**
   * Map xcombs to kernels and generators
   */
  GeneratorMap theGeneratorMap;

  /**
   * The last selected splitting
   */
  Ptr<PowhegSplittingKernel>::tptr theLastSplitting;

  /**
   * True, if verbose.
   */
  bool theVerbose;

  /**
   * True, if the next event should be discarded
   * as the presampling has screwed the parton extractor
   * @TODO Provide PartonExtractor with backup facilities.
   */
  bool theDiscardNext;

  /**
   * Add generators for the newly encountered process
   */
  pair<GeneratorMap::iterator,GeneratorMap::iterator> getGenerators(EventHandler& eh);  

  /**
   * Generate event from given kernel
   */
  Energy generate(pair<Ptr<PowhegSplittingKernel>::ptr,ExponentialGeneratorPtr>&);

  /**
   * Set veto scales, if no radiation generated
   */
  void veto(EventHandler & eh) const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PowhegSplittingGenerator & operator=(const PowhegSplittingGenerator &);

};

}

#endif /* HERWIG_PowhegSplittingGenerator_H */
