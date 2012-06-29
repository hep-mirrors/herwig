// -*- C++ -*-
#ifndef HERWIG_FortranFinder_H
#define HERWIG_FortranFinder_H
//
// This is the declaration of the FortranFinder class.
//

#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Shower/Base/PartnerFinder.h"
#include "FortranFinder.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  typedef of a pair of particle for calculating the evolution scales
 */
typedef pair<tShowerParticlePtr,tShowerParticlePtr> ShowerPPair;

/** \ingroup Shower
 *
 *  The FortranFinder class is responsible for finding the partners and setting the 
 *  initial evolution scales for the original FORTRAN HERWIG shower evolution.
 *
 * @see \ref FortranFinderInterfaces "The interfaces"
 * defined for FortranFinder.
 */
class FortranFinder: public PartnerFinder {

public:

  /**
   * The default constructor.
   */
  inline FortranFinder();

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (Energy(),Energy()) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   */
  //@{
  /**
   *  Calculate the initial evolution scales for two final-state particles
   */
  virtual pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for two initial-state particles
   */
  virtual pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for one initial 
   *  and one final-state particles
   */
  virtual pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
							  const bool isDecayCase);
  //@}

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static NoPIOClassDescription<FortranFinder> initFortranFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FortranFinder & operator=(const FortranFinder &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FortranFinder. */
template <>
struct BaseClassTrait<Herwig::FortranFinder,1> {
  /** Typedef of the first base class of FortranFinder. */
  typedef Herwig::PartnerFinder NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FortranFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FortranFinder>
  : public ClassTraitsBase<Herwig::FortranFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::FortranFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FortranFinder is implemented. It may also include several, space-separated,
   * libraries if the class FortranFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so HwFortranShower.so"; }
};

/** @endcond */

}

#include "FortranFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FortranFinder.tcc"
#endif

#endif /* HERWIG_FortranFinder_H */
