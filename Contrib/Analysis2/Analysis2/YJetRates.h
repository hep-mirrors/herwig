// -*- C++ -*-

// (C) 2007-2009 Simon Plaetzer -- sp@particle.uni-karlsruhe.de

#ifndef Analysis2_YJetRates_H
#define Analysis2_YJetRates_H
//
// This is the declaration of the YJetRates class.
//

#include "Analysis2Base.h"

namespace Analysis2 {

using namespace ThePEG;

/**\ingroup Analysis2
 * 
 * Analysis for jet rates
 *
 * @author Simon Plaetzer
 *
 * @see \ref YJetRatesInterfaces "The interfaces"
 * defined for YJetRates.
 */
class YJetRates: public Analysis2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline YJetRates()
    : _nMax(5), _options(6,""), _yn(), _sumWeights(0.) {}

  /**
   * The destructor.
   */
  virtual ~YJetRates();
  //@}

public:

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Insert options for a given n
   */
  string R (string);

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

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
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}


private:

  /**
   * The maximum n up to which rates should
   * be booked, starting from n=2
   */
  unsigned int _nMax;

  /**
   * A vector of option strings, one for each n
   * The convention is that element 0,1 are ignored,
   * 2 gives options for R2, ...
   */
  vector<string> _options;

  /**
   * Vector of the observable names
   */
  vector<string> _yn;

  /**
   * Sum of all weights booked into
   * jet rates
   */
  double _sumWeights;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<YJetRates> initYJetRates;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  YJetRates & operator=(const YJetRates &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"
#include "ThePEG/Config/Pointers.h"

namespace ThePEG {

ThePEG_DECLARE_POINTERS(Analysis2::YJetRates,YJetRatesPtr);

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of YJetRates. */
template <>
struct BaseClassTrait<Analysis2::YJetRates,1> {
  /** Typedef of the first base class of YJetRates. */
  typedef Analysis2::Analysis2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the YJetRates class and the shared object where it is defined. */
template <>
struct ClassTraits<Analysis2::YJetRates>
  : public ClassTraitsBase<Analysis2::YJetRates> {
  /** Return a platform-independent class name */
  static string className() { return "Analysis2::YJetRates"; }
  /**
   * The name of a file containing the dynamic library where the class
   * YJetRates is implemented. It may also include several, space-separated,
   * libraries if the class YJetRates depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "Analysis2.so"; }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "YJetRates.tcc"
#endif

#endif /* Analysis2_YJetRates_H */
