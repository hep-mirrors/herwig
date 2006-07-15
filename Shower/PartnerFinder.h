// -*- C++ -*-
#ifndef HERWIG_PartnerFinder_H
#define HERWIG_PartnerFinder_H
//
// This is the declaration of the PartnerFinder class.
//

#include "ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "ShowerVariables.h"
#include "PartnerFinder.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  typedef of a pair of particle for calculating the evolution scales
 */
typedef pair<tShowerParticlePtr,tShowerParticlePtr> ShowerPPair;

/** \ingroup Shower
 *
 *  This class is responsible of two related tasks: 
 *  -#   it finds the partners (pairs of ShowerParticle objects)
 *       for each interaction types (QCD,QED,EWK...) and within the
 *       scale range specified by the ShowerVariables object;
 *
 *  -#   for each pair of partners (and interaction therefore)
 *       it sets the initial evolution scales of both of them.
 *  Notice that a given particle has, in general, a different partner
 *  for each different interaction; however, given a partner, its 
 *  initial evolution scale, Q, is purely a kinematical relationship 
 *  between the pair, without dependence on the dynamics (i.e. type of interaction).
 *  Therefore, there must be different methods for finding partners of different
 *  interaction types, but an unique common method to calculate the initial evolving
 *  showering scale.
 *  Notice that:
 *  -    to be more general, one should define an abstract class 
 *       AbsPartnerFinder which has, exactly like the present PartnerFinder
 *       class, the definition of the methods: 
 *       - setQCDInitialEvolutionScales
 *       - setQEDInitialEvolutionScales 
 *       - setEWKInitialEvolutionScales 
 *     which are quite general, whereas the other one is declared pure virtual, 
 *     without definition: 
 *       - virtual ...  calculateInitialEvolutionScales() = 0; 
 *     and then having a concrete PartonFinder class which inherits
 *     from AbsPartnerFinder and provides a definition for this virtual method. 
 *     In fact, it is only in this method that a specific choice of ordering 
 *     variable must be made. Therefore, if we wanted a different choice, 
 *     we could define another class, AlternativePartonFinder, which also 
 *     inherits from AbsPartnerFinder, but provides a different definition 
 *     of the method calculateInitialEvolutionScales.
 *
 *  The following changes to this class need to be considered
 *  - If the method  calculateInitialEvolutionScales
 *    is not used anywhere other than anthoer part of this class,
 *    then it should be moved in the private part of the class.
 *  - It could be easy and straightforward to add to 
 *    setQCDInitialEvolutionScales  
 *    (and, maybe, but less likely, also for the other ones)
 *    either a pointer to the Collision object (which is
 *    already available in the class InsideRangeShowerEvolver
 *    that calls this class Partner Finder) --- in the
 *    case some information from the hard subprocess is
 *    necessary to find colour partners (expecially in the
 *    case of baryon number nonconserving processes in 
 *    Rp violating Susy) --- and/or a boolean flag that
 *    tells whether or not the input particles are
 *    associates to a decay --- and therefore finding
 *    the colour partners, at least for the scale of
 *    the decay, should not involved searches through
 *    the particle history.
 *
 * @see \ref PartnerFinderInterfaces "The interfaces"
 * defined for PartnerFinder.
 */
class PartnerFinder: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline PartnerFinder();
  //@}

public:

  /**
   * Given in input a collection of particles (ShowerParticle objects),
   * each of these methods set the initial evolution scales of those particles, 
   * between the ones given in input, that do not have yet their evolution scale set, 
   * according to a given interaction type (QCD, QED, EWK,...), and within the 
   * constraints specified in the showerVariables object. 
   * The input collection of particles can be either the full collection of 
   * showering particles (kept in the main class ShowerHandler,
   * in the case isDecayCase is false, or simply, in the case isDecayCase 
   * is true, the decaying particle and its decay products.    
   * The methods returns true, unless something wrong (inconsistencies,
   * or undefined values) happens.
   */
  //@{
  /**
   * Set the initial QCD scales
   * @param showerVariables  Pointer to the ShowerVariables
   * @param particles        The particles to be considered
   * @param isDecayCase      Whether or not this is a decay
   */
  bool setQCDInitialEvolutionScales(const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);

  /**
   * Set the initial QED scales
   * @param showerVariables  Pointer to the ShowerVariables
   * @param particles        The particles to be considered
   * @param isDecayCase      Whether or not this is a decay
   */
  bool setQEDInitialEvolutionScales(const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);

  /**
   * Set the initial electroweak scales
   * @param showerVariables  Pointer to the ShowerVariables
   * @param particles        The particles to be considered
   * @param isDecayCase      Whether or not this is a decay
   */
  bool setEWKInitialEvolutionScales(const ShowerParticleVector &particles,
                                    const bool isDecayCase = false);
  //@}

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (Energy(),Energy()) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   */
  //@{
  /**
   *  General method to calculate the initial evolution scales
   */
  pair<Energy,Energy> calculateInitialEvolutionScales(const ShowerPPair &,
                                           const bool isDecayCase = false);

  /**
   *  Calculate the initial evolution scales for two final-state particles
   */
  pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for two initial-state particles
   */
  pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &);

  /**
   *  Calculate the initial evolution scales for one initial 
   *  and one final-state particles
   */
  pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
                                       const bool isDecayCase = false);
  //@}

  /**
   *  Set the shower variables, only used by Evolver in doinit
   */
  inline void setShowerVariables(ShowerVarsPtr);

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
  static ClassDescription<PartnerFinder> initPartnerFinder;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PartnerFinder & operator=(const PartnerFinder &);

private:

  /**
   *  Approach to use for setting the colour partners
   */
  int _approach;

  /**
   *  Pointer to the ShowerVariables object
   */
  ShowerVarsPtr _showerVariables;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of PartnerFinder. */
template <>
struct BaseClassTrait<Herwig::PartnerFinder,1> {
  /** Typedef of the first base class of PartnerFinder. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the PartnerFinder class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::PartnerFinder>
  : public ClassTraitsBase<Herwig::PartnerFinder> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::PartnerFinder"; }
  /**
   * The name of a file containing the dynamic library where the class
   * PartnerFinder is implemented. It may also include several, space-separated,
   * libraries if the class PartnerFinder depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "PartnerFinder.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "PartnerFinder.tcc"
#endif

#endif /* HERWIG_PartnerFinder_H */
