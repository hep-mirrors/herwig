// -*- C++ -*-
#ifndef HERWIG_SimpleILCAnalysis_H
#define HERWIG_SimpleILCAnalysis_H
//
// This is the declaration of the SimpleILCAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the SimpleILCAnalysis class.
 *
 * @see \ref SimpleILCAnalysisInterfaces "The interfaces"
 * defined for SimpleILCAnalysis.
 */
class SimpleILCAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SimpleILCAnalysis();

 public:
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  //static void Init();

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

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  //  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  //  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  // virtual void analyze(tPPtr particle);
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  //  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  // void persistentInput(PersistentIStream & is, int version);
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<SimpleILCAnalysis> initSimpleILCAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleILCAnalysis & operator=(const SimpleILCAnalysis &);

private:
 
   /**
   *  The weight for the event
   */
  double eventweight_;
 /**
   * Mass of the top pair
   */
  Histogram _Mt;
  Histogram _Mtb;

  /**
   *  Momenta of the b, bbar (energy pbe,pbbe, transverse mmt 
    *ptb ptbb, long. mmt pzb, pzbb and ratios to the energies 
   * ending in *r
   */
  Histogram _pbe; 
  Histogram _pbbe;
  Histogram _ptb;
  Histogram _ptbb;
  Histogram _ptbr;
  Histogram _ptbbr;
  Histogram _pzb;
  Histogram _pzbb;
  Histogram _pzbr;
  Histogram _pzbbr;
    /**
   * Rapidity of b, bbar
   */
  Histogram _rapb;
  Histogram _rapbb;
   /**
   * Correlation angles
   */
  Histogram _acostem;  // btw Decay electron and top
  Histogram _acostee;  //  btw Decay charged leptons
  Histogram _acostebe; // btw Decay electron and incoming e-
  

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SimpleILCAnalysis. */
template <>
struct BaseClassTrait<Herwig::SimpleILCAnalysis,1> {
  /** Typedef of the first base class of SimpleILCAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SimpleILCAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SimpleILCAnalysis>
  : public ClassTraitsBase<Herwig::SimpleILCAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SimpleILCAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SimpleILCAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "SimpleILCAnalysis.so"; }
};

/** @endcond */

}

//#include "SimpleILCAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SimpleILCAnalysis.tcc"
#endif

#endif /* HERWIG_SimpleILCAnalysis_H */
