// -*- C++ -*-
#ifndef HERWIG_MySimpleAnalysis_H
#define HERWIG_MySimpleAnalysis_H
//
// This is the declaration of the MySimpleAnalysis class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "Herwig++/Utilities/Histogram.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MySimpleAnalysis class.
 *
 * @see \ref MySimpleAnalysisInterfaces "The interfaces"
 * defined for MySimpleAnalysis.
 */
class MySimpleAnalysis: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MySimpleAnalysis();

  /**
   * The destructor.
   */
  //virtual void analyze(tEventPtr event, long ieve, int loop, int state);
  //virtual ~MySimpleAnalysis();
  //@}

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

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  // virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  //virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  //virtual void analyze(tPPtr particle);
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  //void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  //void persistentInput(PersistentIStream & is, int version);
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
  static NoPIOClassDescription<MySimpleAnalysis> initMySimpleAnalysis;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MySimpleAnalysis & operator=(const MySimpleAnalysis &);

private:
  /**
   *   \f$p_T\f$ of the Z boson
   */
  vector<Histogram> _ptZ;
  vector<Histogram> _ptWp;
  vector<Histogram> _ptWm;
 
 /**
   *   \f$p_T\f$ of the W/Z boson for data comparism
   */
  Histogram _ptZ4,_ptZ5,_ptZ6,_ptZ7,_ptZ8,_ptZ9,_ptZ10;
  Histogram _ptW4,_ptW5,_ptW6,_ptW7,_ptW8,_ptW9,_ptW10;
  
  /**
   * Mass of the Z boson
   */
  Histogram _mZ;
  Histogram _mWp;
  Histogram _mWm;
  /**
   *  The weight for the event
   */
  double eventweight_;
  /**
   *  Rapidity of Z
   */
  Histogram _rapZ;
  Histogram _rapWp;
  Histogram _rapWm;

  /**
   *  Azimuth of Z
   */
  Histogram _phiZ;
  Histogram _phiWp;
  Histogram _phiWm;
/**
   *  Rapidity of leptons
   */
  Histogram _rapem;
  Histogram _rapep;
  Histogram _rapnu;
  Histogram _rapanu;
  /**
   *  Azimuth of leptons
   */
  Histogram _phiem;
  Histogram _phiep;
  Histogram _phinu;
  Histogram _phianu;
  
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MySimpleAnalysis. */
template <>
struct BaseClassTrait<Herwig::MySimpleAnalysis,1> {
  /** Typedef of the first base class of MySimpleAnalysis. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MySimpleAnalysis class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MySimpleAnalysis>
  : public ClassTraitsBase<Herwig::MySimpleAnalysis> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MySimpleAnalysis"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MySimpleAnalysis class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "MySimpleAnalysis.so"; }
};

/** @endcond */

}

//#include "MySimpleAnalysis.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MySimpleAnalysis.tcc"
#endif

#endif /* HERWIG_MySimpleAnalysis_H */
