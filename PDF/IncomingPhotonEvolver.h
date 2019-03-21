// -*- C++ -*-
#ifndef HERWIG_IncomingPhotonEvolver_H
#define HERWIG_IncomingPhotonEvolver_H
//
// This is the declaration of the IncomingPhotonEvolver class.
//

#include "ThePEG/Handlers/StepHandler.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The IncomingPhotonEvolver class performs the backward evolution
 * of a photon in a partonic process to a quark or antiquark so that
 * the event can be showered.
 *
 * @see \ref IncomingPhotonEvolverInterfaces "The interfaces"
 * defined for IncomingPhotonEvolver.
 */
class IncomingPhotonEvolver: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  IncomingPhotonEvolver();
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

  /**
   * Compute boost parameter along z axis to get (Ep, any perp, qp)
   * from (E, same perp, q).
   */
  inline double getBeta(const double E, const double q, 
			const double Ep, const double qp) const
  {return (q*E-qp*Ep)/(sqr(qp)+sqr(E));}

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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<IncomingPhotonEvolver> initIncomingPhotonEvolver;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IncomingPhotonEvolver & operator=(const IncomingPhotonEvolver &) = delete;

private:
  /**
   * PDF set to use. Overrides the one that is associated with the beam particle.
   */
  PDFPtr PDF_;

  /**
   *  The maximum value of the PDF for the sample
   */
  double PDFMax_;

  /**
   *  The power for the sampling of the PDF
   */
  double PDFPower_;

  /**
   *  The minimum starting scale for the evolution
   */
  Energy minpT_;

  /**
   *  The minimum space-like virtuality of the photon
   */
  Energy minVirtuality_;

  /**
   *  Maximum number of attempts to generate the scale of the
   *  branching
   */
  unsigned int vetoTries_;

  /**
   *  Maximum number of attempts to regenerate the virtuality
   */
  unsigned int virtualityTries_;

  /**
   *   Photon ParticleData object
   */
  tcPDPtr photon_;

  /**
   *  Partons to backward evolve to
   */
  vector<tcPDPtr> partons_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of IncomingPhotonEvolver. */
template <>
struct BaseClassTrait<Herwig::IncomingPhotonEvolver,1> {
  /** Typedef of the first base class of IncomingPhotonEvolver. */
  typedef StepHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the IncomingPhotonEvolver class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::IncomingPhotonEvolver>
  : public ClassTraitsBase<Herwig::IncomingPhotonEvolver> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::IncomingPhotonEvolver"; }
  /**
   * The name of a file containing the dynamic library where the class
   * IncomingPhotonEvolver is implemented. It may also include several, space-separated,
   * libraries if the class IncomingPhotonEvolver depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwIncomingPhotonEvolver.so"; }
};

/** @endcond */

}

#endif /* HERWIG_IncomingPhotonEvolver_H */
