// -*- C++ -*-
#ifndef Herwig_BranchingRatioReweighter_H
#define Herwig_BranchingRatioReweighter_H
//
// This is the declaration of the BranchingRatioReweighter class.
//

#include "ThePEG/Handlers/StepHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Decay
 * The BranchingRatioReweighter class is designed to reweight events
 * where some decay modes of a particle, or many particles, have been
 * switched off in order to improve the statistics.
 *
 * @see \ref BranchingRatioReweighterInterfaces "The interfaces"
 * defined for BranchingRatioReweighter.
 */
class BranchingRatioReweighter: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  BranchingRatioReweighter();

  /**
   * The destructor.
   */
  virtual ~BranchingRatioReweighter();
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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BranchingRatioReweighter & operator=(const BranchingRatioReweighter &);

};

}

#endif /* Herwig_BranchingRatioReweighter_H */
