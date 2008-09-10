// -*- C++ -*-
//
// MECorrectionBase.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MECorrectionBase_H
#define HERWIG_MECorrectionBase_H
//
// This is the declaration of the MECorrectionBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "ShowerProgenitor.h"
#include "ShowerTree.fh"
#include "Evolver.fh"
#include "MECorrectionBase.fh"

namespace Herwig {

using namespace ThePEG;

/**\ingroup Shower
 *
 * This is the base class for all matrix element corrections
 * used inside parton showers.
 *
 * @see \ref MECorrectionBaseInterfaces "The interfaces"
 * defined for MECorrectionBase.
 */
class MECorrectionBase: public Interfaced {

/**
 *  The Evolver is a friend to set some variables at initialisation
 */
friend class Evolver;

public:

  /**
   *  Virtual members which should be implemented in classes inheriting
   * from this one in order to apply the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   * @param tree The shower tree currently being showered
   * @param initial The initial-state radiation enhancement factor
   * @param final   The final-state radiation enhancement factor
   * @param evolver Pointer to the Evolver.
   */
  virtual bool canHandle(ShowerTreePtr tree,double & initial,
			 double & final,EvolverPtr evolver)=0;

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr)=0;

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br)=0;
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
   *  Set/Get methods for the pointer to the evolver
   */
  //@{
  /**
   *  Set the evolver
   */
  void evolver(tEvolverPtr);

  /**
   *  Get the evolver
   */
  tEvolverPtr evolver() const;
  //@}

  /**
   *  Access to the coupling
   */
  ShowerAlphaPtr coupling() {return _alpha;}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractClassDescription<MECorrectionBase> initMECorrectionBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MECorrectionBase & operator=(const MECorrectionBase &);

private:

  /**
   *  Pointer to the coupling
   */
  ShowerAlphaPtr _alpha;

  /**
   *  Pointer to the Evolver object
   */
  tEvolverPtr _evolver;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MECorrectionBase. */
template <>
struct BaseClassTrait<Herwig::MECorrectionBase,1> {
  /** Typedef of the first base class of MECorrectionBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MECorrectionBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MECorrectionBase>
  : public ClassTraitsBase<Herwig::MECorrectionBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MECorrectionBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MECorrectionBase is implemented. It may also include several, space-separated,
   * libraries if the class MECorrectionBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MECorrectionBase_H */
