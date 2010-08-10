// -*- C++ -*-
#ifndef HERWIG_HiddenPartnerFinder_H
#define HERWIG_HiddenPartnerFinder_H
//
// This is the declaration of the HiddenPartnerFinder class.
//

#include "Herwig++/Shower/Default/QTildeFinder.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the HiddenPartnerFinder class.
 *
 * @see \ref HiddenPartnerFinderInterfaces "The interfaces"
 * defined for HiddenPartnerFinder.
 */
class HiddenPartnerFinder: public QTildeFinder {

public:

  /**
   * Given in input a collection of particles (ShowerParticle objects),
   * each of these methods set the initial evolution scales of those particles, 
   * between the ones given in input, that do not have yet their
   * evolution scale set. 
   * The input collection of particles can be either the full collection of 
   * showering particles (kept in the main class ShowerHandler,
   * in the case isDecayCase is false, or simply, in the case isDecayCase 
   * is true, the decaying particle and its decay products.    
   * The methods returns true, unless something wrong (inconsistencies,
   * or undefined values) happens.
   *
   * These methods are virtual but in most cases inheriting classes should not
   * need to overide them as they simply find the relevant partner and call
   * one of the calculateScale members to calculate the scale.
   */
  //@{
  /**
   * Set the initial scales
   * @param particles        The particles to be considered
   * @param isDecayCase      Whether or not this is a decay
   * @param setPartners Whether to set the colour partners or just the scales
   */
  virtual bool setInitialEvolutionScales(const ShowerParticleVector &particles,
					 const bool isDecayCase,
					 const bool setPartners=true);
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
  static ClassDescription<HiddenPartnerFinder> initHiddenPartnerFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HiddenPartnerFinder & operator=(const HiddenPartnerFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HiddenPartnerFinder. */
template <>
struct BaseClassTrait<Herwig::HiddenPartnerFinder,1> {
  /** Typedef of the first base class of HiddenPartnerFinder. */
  typedef Herwig::QTildeFinder NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HiddenPartnerFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::HiddenPartnerFinder>
  : public ClassTraitsBase<Herwig::HiddenPartnerFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::HiddenPartnerFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HiddenPartnerFinder is implemented. It may also include several, space-separated,
   * libraries if the class HiddenPartnerFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HiddenValleyModel.so"; }
};

/** @endcond */

}

#endif /* HERWIG_HiddenPartnerFinder_H */
