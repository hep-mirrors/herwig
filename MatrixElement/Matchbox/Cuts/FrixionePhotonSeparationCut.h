// -*- C++ -*-
//
// FrixionePhotonSeparationCut.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_FrixionePhotonSeparationCut_H
#define Herwig_FrixionePhotonSeparationCut_H
//
// This is the declaration of the FrixionePhotonSeparationCut class.
//

#include "ThePEG/Cuts/MultiCutBase.h"
#include "ThePEG/PDT/MatcherBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Michael Rauch
 *
 * \brief This class implements a cut on legoplot and rapidity separation
 *
 * @see \ref FrixionePhotonSeparationCutInterfaces "The interfaces"
 * defined for FrixionePhotonSeparationCut.
 */
class FrixionePhotonSeparationCut: public MultiCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FrixionePhotonSeparationCut() :  theDeltaZero(0.0), theExponentn(1.0), theEfficiency(0.0), theCutType(1) {}
  //@}

public:

  /** @name Overridden virtual functions defined in the base class. */
  //@{
  /**
   * Return true if a pair of particles with type \a pitype and \a
   * pjtype and momenta \a pi and \a pj respectively passes the
   * cuts. \a inci and \a inj indicates if the corresponding particles
   * are incoming.
   */
  virtual bool passCuts(tcCutsPtr, const tcPDVector & ptype,
	                const vector<LorentzMomentum> & p) const;
  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

  /**
   * Return the matcher for particles to isolate on.
   */
  Ptr<MatcherBase>::tptr matcher() const { return theMatcher; }

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

private:

  /**
   * The maximal legoplot separation where partons are included in the criterium
   */
  double theDeltaZero;

  /**
   * The exponent n of the algorithm
   */
  double theExponentn;

  /**
   * The efficiency epsilon of the algorithm
   */
  double theEfficiency;

  /**
   * The cut type of the algorithm
   */
  int theCutType;

  /**
   * A matcher for particles to isolate on.
   */
  Ptr<MatcherBase>::ptr theMatcher;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<FrixionePhotonSeparationCut> initFrixionePhotonSeparationCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FrixionePhotonSeparationCut & operator=(const FrixionePhotonSeparationCut &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FrixionePhotonSeparationCut. */
template <>
struct BaseClassTrait<Herwig::FrixionePhotonSeparationCut,1> {
  /** Typedef of the first base class of FrixionePhotonSeparationCut. */
  typedef MultiCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FrixionePhotonSeparationCut class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FrixionePhotonSeparationCut>
  : public ClassTraitsBase<Herwig::FrixionePhotonSeparationCut> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::FrixionePhotonSeparationCut"; }
  /** Return the name of the shared library be loaded to get
   *  access to the FrixionePhotonSeparationCut class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwMatchboxCuts.so"; }
};

/** @endcond */

}

#endif /* Herwig_FrixionePhotonSeparationCut_H */
