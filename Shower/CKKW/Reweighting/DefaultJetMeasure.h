// -*- C++ -*-
//
// DefaultJetMeasure.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DefaultJetMeasure_H
#define HERWIG_DefaultJetMeasure_H
//
// This is the declaration of the DefaultJetMeasure class.
//

#include "JetMeasure.h"
#include "DefaultJetMeasure.fh"

#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**\ingroup CKKW
 *
 * DefaultJetMeasure is a jet resolution for the standard
 * CKKW approaches together with the standard shower implementation.
 *
 * It provides methods to veto shower emissions, returns phase
 * space boundaries as well as Jacobians for Sudakov exponent
 * integration.
 *
 *@author Simon Plaetzer
 *
 * @see \ref DefaultJetMeasureInterfaces "The interfaces"
 * defined for DefaultJetMeasure.
 */
class DefaultJetMeasure: public JetMeasure {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~DefaultJetMeasure();
  //@}

protected:

  /**
   * Construct giving the actual type of shower.
   * num should be two for 1->2 showers, and 3
   * for showers with a Sudakov decomposition,
   * i.e. when a spectator is to be considered
   * in evaluating the kinematics.
   */
  explicit inline DefaultJetMeasure (unsigned int num);

public:

  /**
   * Set the sudakov basis to be used
   */
  inline void sudakovBasis (const Lorentz5Momentum& p, const Lorentz5Momentum& n);
  
  /**
   * Get the sudakov basis to be used
   */
  inline pair<Lorentz5Momentum,Lorentz5Momentum> sudakovBasis () const;

  /**
   * Return true, if the jet resolution can
   * handle the given configuration.
   */
  virtual bool canHandle (const IdList&, bool initial = false) = 0;

  /**
   * Return true, if the branching configuration
   * together with the scale and momentum fraction
   * is resolvable. This is used for integrating
   * Sudakov exponents.
   */
  virtual bool resolvable (const IdList&, Energy2, Energy2, double, bool initial = false) = 0;

  /**
   * Return the Jacobian for going from the parton shower
   * variables to the clustering variables as a function
   * of the clustering variables.
   */
  virtual double showerJacobian (const IdList&, Energy2, double, bool initial = false) = 0;

  /**
   * Return the shower scale and momentum fraction
   * as a function of the clustering variables.
   */
  virtual pair<Energy2,double> invertClustering (const IdList&, Energy2, double, bool initial = false) = 0;


  /**
   * Return true, if the given branching is resolvable.
   * This is used for vetoing shower emissions. An
   * optional resolution scale may be supplied. If this
   * is set to zero, the predefined resolution is to be used.
   */
  virtual bool resolvable (tcShowerParticlePtr,
			   const Branching&,
			   bool initial = false,
			   Energy2 resolution = 0.*GeV2) = 0;

  /**
   * Return the bounds on the clustering momentum
   * fractions as a function of scale and the hard
   * scale.
   */
  virtual pair<double,double> zLimits (Energy2 q, Energy2 Q, const IdList&, bool initial = false) = 0;


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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The sudakov basis to be used
   */
  pair<Lorentz5Momentum,Lorentz5Momentum> _sudakovBasis;

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DefaultJetMeasure> initDefaultJetMeasure;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DefaultJetMeasure & operator=(const DefaultJetMeasure &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DefaultJetMeasure. */
template <>
struct BaseClassTrait<Herwig::DefaultJetMeasure,1> {
  /** Typedef of the first base class of DefaultJetMeasure. */
  typedef Herwig::JetMeasure NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DefaultJetMeasure class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DefaultJetMeasure>
  : public ClassTraitsBase<Herwig::DefaultJetMeasure> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DefaultJetMeasure"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DefaultJetMeasure is implemented. It may also include several, space-separated,
   * libraries if the class DefaultJetMeasure depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "DefaultJetMeasure.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DefaultJetMeasure.tcc"
#endif

#endif /* HERWIG_DefaultJetMeasure_H */
