// -*- C++ -*-
#ifndef Herwig_SpinHadronizer_H
#define Herwig_SpinHadronizer_H
//
// This is the declaration of the SpinHadronizer class.
//

#include "ThePEG/Handlers/StepHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The SpinHadronizer class is designed to be used as a post-hadronization handler to
 * give a simple model of spin transfer between the perturbative and non-perturbative
 * stages.
 *
 * @see \ref SpinHadronizerInterfaces "The interfaces"
 * defined for SpinHadronizer.
 */
class SpinHadronizer: public StepHandler {

public:

  /**
   * The default constructor.
   */
  SpinHadronizer() : omegaHalf_(2./3.), omegaThreeHalf_(0.2),
        minFlav_(3), maxFlav_(5), debug_(false), qPol_(6,make_pair(0.,0.))
  {}

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

  /**
   *  Functions to calculate the spins
   */
  //@{

  /**
   *  Calculate the spin of a baryon
   */
  void baryonSpin(tPPtr baryon);

  /**
   *  Calculate the spin of a meson
   */
  void mesonSpin(tPPtr meson);

  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

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
  SpinHadronizer & operator=(const SpinHadronizer &);

private:

  /**
   *  Parameters
   */
  //@{
  /**
   *  Falk-Peskin \f$\omega_\frac12\f$ parameter
   */
  double omegaHalf_;

  /**
   *  Falk-Peskin \f$\omega_\frac32\f$ parameter
   */
  double omegaThreeHalf_;

  /**
   *  Minimum quark flavour
   */
  unsigned int minFlav_;

  /**
   *  Maximum quark flavour
   */
  unsigned int maxFlav_;

  /**
   *  Print out debugging info
   */
  bool debug_;

  /**
   *  Polarization of the quarks
   */
  vector<pair<double,double> > qPol_;
  //@}

};

}

#endif /* Herwig_SpinHadronizer_H */
