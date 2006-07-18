// -*- C++ -*-
#ifndef HERWIG_MultiplicityCount_H
#define HERWIG_MultiplicityCount_H
//
// This is the declaration of the MultiplicityCount class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/ParticleData.h"
#include "MultiplicityCount.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  Enumeration for species of particle
 */
enum ParticleSpecies 
{
  lightMeson=0,strangeMeson,lightBaryon,other
};

/**
 *  Struct for the multiplcity data
 */
struct MultiplicityInfo
{
  /**
   *  Default constructor
   * @param mult  The observed multiplcity.
   * @param error The error on the observed multiplicity
   * @param type  The type of particle
   */
  inline MultiplicityInfo(double mult=0.,double error=0.,
			  ParticleSpecies type=other);

  /**
   *  The observed multiplicity
   */
  double mult;

  /**
   *  The error on the observed multiplicity
   */
  double error;

  /**
   *  The type of particle
   */
  ParticleSpecies type;

  /**
   *  Number of particles of this type
   */
  long actualCount;

  /**
   *  Sum of squares of number per event for error
   */
  double sumofsquares;

  /**
   *  The average number per event
   * @param N The number of events
   */
  double mean(double N);

  /**
   *  The error on the average number per event
   * @param N The number of events 
   */
  double stderror(double N);

  /**
   *  Is the result more than \f$3\sigma\f$ from the experimental result
   * @param N The number of events
   */
  bool serious(double N);
};

/**
 * The MultiplicityCount class is designed to count particle multiplicities and
 * compare them to LEP data.
 *
 * @see \ref MultiplicityCountInterfaces "The interfaces"
 * defined for MultiplicityCount.
 */
class MultiplicityCount: public AnalysisHandler {

public:

  /**
   * The default constructor.
   */
  MultiplicityCount();

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
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
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
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<MultiplicityCount> initMultiplicityCount;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiplicityCount & operator=(const MultiplicityCount &);

private:

  /**
   *  The PDG codes of the particles
   */
  vector<long> _particlecodes;

  /**
   * The multiplcity
   */
  vector<double> _multiplicity;

  /**
   * The error
   */
  vector<double> _error;

  /**
   * Species of particle
   */
  vector<unsigned int> _species;

  /**
   *  Map of PDG codes to multiplicity info
   */
  map<long,MultiplicityInfo> _data;

  /**
   *  Map of number of final-state particles to PDG code
   */
  map<long,long> _finalstatecount;

  /**
   *  Particles in hard process
   */
  map<long,long> _collisioncount;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MultiplicityCount. */
template <>
struct BaseClassTrait<Herwig::MultiplicityCount,1> {
  /** Typedef of the first base class of MultiplicityCount. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MultiplicityCount class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MultiplicityCount>
  : public ClassTraitsBase<Herwig::MultiplicityCount> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::MultiplicityCount"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the MultiplicityCount class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwAnalysis.so"; }
};

/** @endcond */

}

#include "MultiplicityCount.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "MultiplicityCount.tcc"
#endif

#endif /* HERWIG_MultiplicityCount_H */
