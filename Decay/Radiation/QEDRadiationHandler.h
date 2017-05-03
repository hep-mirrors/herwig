// -*- C++ -*-
//
// QEDRadiationHandler.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_QEDRadiationHandler_H
#define HERWIG_QEDRadiationHandler_H
//
// This is the declaration of the QEDRadiationHandler class.
//

#include "ThePEG/Handlers/StepHandler.h"
#include "DecayRadiationGenerator.fh"
#include "QEDRadiationHandler.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The QEDRadiationHandler class is designed so that the approach
 * for the generation of QED radiation in decays can be used to generate
 * QED radiation from the decay products of \f$s\f$-channel resonances
 * produced in the hard process where the decay products of the resonance
 * have been generated as part of the hard process.
 *
 * It is designed to be used as a PostSubProcessHandler.
 *
 * @see \ref QEDRadiationHandlerInterfaces "The interfaces"
 * defined for QEDRadiationHandler.
 */
class QEDRadiationHandler: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QEDRadiationHandler();
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
  static ClassDescription<QEDRadiationHandler> initQEDRadiationHandler;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QEDRadiationHandler & operator=(const QEDRadiationHandler &);

private:

  /**
   *  Pointer to the object responsible for generating the radiation in
   *  the decays
   */
  DecayRadiationGeneratorPtr _generator;

  /**
   *  List of the PDG codes of the decaying particles which should be considered
   */
  vector<long> _decayingParticles;

  /**
   *  List of the PDG codes of the decay products which should be considered
   */
  vector<long> _decayProducts;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of QEDRadiationHandler. */
template <>
struct BaseClassTrait<Herwig::QEDRadiationHandler,1> {
  /** Typedef of the first base class of QEDRadiationHandler. */
  typedef StepHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the QEDRadiationHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::QEDRadiationHandler>
  : public ClassTraitsBase<Herwig::QEDRadiationHandler> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::QEDRadiationHandler"; }
};

/** @endcond */

}

#endif /* HERWIG_QEDRadiationHandler_H */
