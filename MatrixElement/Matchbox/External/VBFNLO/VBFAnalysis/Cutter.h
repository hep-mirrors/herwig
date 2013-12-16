// -*- C++ -*-
#ifndef THEPEG_Cutter_H
#define THEPEG_Cutter_H
//
// This is the declaration of the Cutter class.
//

#include "ThePEG/Handlers/StepHandler.h"
#include "NLOAnalysis.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "BasicHistogramCollection.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the Cutter class.
 *
 * @see \ref CutterInterfaces "The interfaces"
 * defined for Cutter.
 */
class Cutter: public StepHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Cutter();

  /**
   * The destructor.
   */
  virtual ~Cutter();
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
		      const Hint & hint) throw(Veto, Stop, Exception);
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<Cutter> initCutter;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Cutter & operator=(const Cutter &);

  
  //AnaPtr Excessevanahandler() {return theexcessevanahandler;}
  //  tcAnaPtr FJAnalysis() {return theFJAnalysis;}
  //  tcAnaPtr Excessevanahandler() {return theexcessevanahandler;}

  Ptr<NLOAnalysis>::ptr theNLOAnalysis;
  Ptr<NLOAnalysis>::ptr theExcessEventsAnalysis;
  //  Ptr<NLOAnalysis>::ptr theexcessevanahandler;

public:
  
  /**
   * Check if the particle is allowed for jet reconstruction
   */
  bool allowedParticle (const Particle &p);
  
  /**
   * Push particles into a PseudoJet vector with the option of sorting out
   */
  vector<fastjet::PseudoJet> recombinables(const ParticleVector& p, bool sort_out = true);

  /**
   * Push particles into a PseudoJet vector with the option of sorting out
   */
  vector<fastjet::PseudoJet> recombinables(const tPVector& p, bool sort_out = true);

  /**
   * Recombine to jets using fastjet methods
   */
  vector<fastjet::PseudoJet> recombine(const vector<fastjet::PseudoJet>& p);

  /**
   * Return jets within the detector range
   */
  vector<fastjet::PseudoJet> getInRange(const vector<fastjet::PseudoJet>& j);

  /*
  double Rparam() const;
  double eta_max() const;//=5.0;
  double Rjj_min() const;//=0.6;;
  double pt_min() const;// = 30.0;
  */

  double theRparam;
  double theEtaDetector;
  double theJetDefPt;


  bool passCuts(vector<fastjet::PseudoJet>);
  //  bool passNjetCut(vector<fastjet::PseudoJet>);
  //bool passDelYhjCut(vector<fastjet::PseudoJet>, fastjet::PseudoJet);
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Cutter. */
template <>
struct BaseClassTrait<arnold::Cutter,1> {
  /** Typedef of the first base class of Cutter. */
  typedef StepHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Cutter class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::Cutter>
  : public ClassTraitsBase<Herwig::Cutter> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::Cutter"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Cutter is implemented. It may also include several, space-separated,
   * libraries if the class Cutter depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "VBFAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_Cutter_H */
