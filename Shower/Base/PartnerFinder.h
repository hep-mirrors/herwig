// -*- C++ -*-
//
// PartnerFinder.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_PartnerFinder_H
#define HERWIG_PartnerFinder_H
//
// This is the declaration of the PartnerFinder class.
//
#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Evolver.fh"
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
 *  -  it finds the partners 
 *  -  for each pair of partners (and interaction therefore)
 *      it sets the initial evolution scales of both of them.
 *
 *  In general the finding of the partners is performed by this class but
 *  the calculation of the initial evolution scales should be implemented
 *  for different shower evolution models in classes inheriting from this one.
 *  Notice that a given particle has, in general, a different partner
 *  for each different interaction; however, given a partner, its 
 *  initial evolution scale, Q, is purely a kinematical relationship 
 *  between the pair, without dependence on the dynamics (i.e. type of interaction).
 *
 * @see \ref PartnerFinderInterfaces "The interfaces"
 * defined for PartnerFinder.
 */
class PartnerFinder: public Interfaced {

public:

  /**
   * The default constructor.
   */
  PartnerFinder() : _approach(0) {}

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

  /**
   * Given a pair of particles, supposedly partners w.r.t. an interaction,
   * this method returns their initial evolution scales as a pair.
   * If something wrong happens, it returns the null (ZERO,ZERO) pair. 
   * This method is used by the above setXXXInitialEvolutionScales 
   * methods.
   * These methods must be overiden in inheriting classes
   */
  //@{
  /**
   *  General method to calculate the initial evolution scales
   */
  virtual pair<Energy,Energy> calculateInitialEvolutionScales(const ShowerPPair &,
							      const bool isDecayCase);

  /**
   *  Calculate the initial evolution scales for two final-state particles
   */
  virtual pair<Energy,Energy> calculateFinalFinalScales(const ShowerPPair &)=0;

  /**
   *  Calculate the initial evolution scales for two initial-state particles
   */
  virtual pair<Energy,Energy> calculateInitialInitialScales(const ShowerPPair &)=0;

  /**
   *  Calculate the initial evolution scales for one initial 
   *  and one final-state particles
   */
  virtual pair<Energy,Energy> calculateInitialFinalScales(const ShowerPPair &,
							  const bool isDecayCase)=0;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractClassDescription<PartnerFinder> initPartnerFinder;

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
  static string className() { return "Herwig::PartnerFinder"; }
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

#endif /* HERWIG_PartnerFinder_H */
