// -*- C++ -*-
//
// DipoleSplittingReweight.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleSplittingReweight_H
#define HERWIG_DipoleSplittingReweight_H
//
// This is the declaration of the DipoleSplittingReweight class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleSplittingInfo.h"

#include "Herwig/Shower/Dipole/DipoleShowerHandler.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleSplittingReweight is used by the dipole shower
 * to reweight splittings from a given dipole splitting kernel.
 *
 * @see \ref DipoleSplittingReweightInterfaces "The interfaces"
 * defined for DipoleSplittingReweight.
 */
class DipoleSplittingReweight: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleSplittingReweight();

  /**
   * The destructor.
   */
  virtual ~DipoleSplittingReweight();
  //@}

public:

  /**
   * Return true, if the reweighting should be applied to the first
   * interaction
   */
  virtual bool firstInteraction() const { return true; }

  /**
   * Return true, if the reweighting should be applied to the secondary
   * interactions
   */
  virtual bool secondaryInteractions() const { return false; }

  /**
   * Update the pointer to the currently active dipole shower handler object.
   */
  void updateCurrentHandler();

  /**
   * Return the pointer to the currently active dipole shower handler object.
   */
  Ptr<DipoleShowerHandler>::tptr currentHandler() const;

  /**
   * Return the reweighting factor for the given splitting type.
   */
  virtual double evaluate(const DipoleSplittingInfo&) const = 0;

  /**
   * Return an enhancement hint for the sampling of the un-reweighted
   * splitting kernel
   */
  virtual double hint(const DipoleSplittingInfo&) const {
    return 1.;
  }

  /**
   * Return true, if the reweight can be entirely absorbed into the hint. A
   * possible detuning will be switched off.
   */
  virtual bool hintOnly(const DipoleSplittingInfo&) const {
    return false;
  }

  /**
   * Set the factor in front of enhance used by the veto algorithm.
   */
  virtual void reweightFactor(const double) {
    return;
  }

  /**
   * Scaling factor for negative reweights.
   */
  virtual void negativeScaling(const double) {
    return;
  }

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
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static AbstractClassDescription<DipoleSplittingReweight> initDipoleSplittingReweight;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleSplittingReweight & operator=(const DipoleSplittingReweight &) = delete;

  /**
   * A pointer to the currently active dipole shower handler object.
   */
  Ptr<DipoleShowerHandler>::tptr theCurrentHandler;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleSplittingReweight. */
template <>
struct BaseClassTrait<Herwig::DipoleSplittingReweight,1> {
  /** Typedef of the first base class of DipoleSplittingReweight. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleSplittingReweight class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleSplittingReweight>
  : public ClassTraitsBase<Herwig::DipoleSplittingReweight> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleSplittingReweight"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleSplittingReweight is implemented. It may also include several, space-separated,
   * libraries if the class DipoleSplittingReweight depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleSplittingReweight_H */
