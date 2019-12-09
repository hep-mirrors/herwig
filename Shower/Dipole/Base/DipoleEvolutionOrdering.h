// -*- C++ -*-
//
// DipoleEvolutionOrdering.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleEvolutionOrdering_H
#define HERWIG_DipoleEvolutionOrdering_H
//
// This is the declaration of the DipoleEvolutionOrdering class.
//

#include "ThePEG/Handlers/HandlerBase.h"

#include "Herwig/Shower/Dipole/Base/Dipole.h"
#include "Herwig/Shower/Dipole/Base/DipoleChain.h"
#include "Herwig/Shower/Dipole/Kernels/DipoleSplittingKernel.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleEvolutionOrdering defines a particular evolution
 * algortihm for the dipole shower.
 *
 */
class DipoleEvolutionOrdering: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleEvolutionOrdering();

  /**
   * The destructor.
   */
  virtual ~DipoleEvolutionOrdering();
  //@}

public:

  /**
   * Return true, if this ordering requests
   * an independent dipole evolution.
   */
  virtual bool independentDipoles() const { return false; }

  /**
   * For the given dipole and splitting kernel return
   * the hard scale.
   */
  virtual Energy hardScale(tPPtr emitter, tPPtr spectator,
			   double emitterX, double spectatorX,
			   const DipoleSplittingKernel&,
			   const DipoleIndex&) const = 0;

  /**
   * For the given performed splitting, dipole chain
   * and dipoles originating from the splitting, set the next
   * scale.
   */
  virtual void setEvolutionScale(Energy scale,
				 const DipoleSplittingInfo&,
				 DipoleChain&,
				 pair<list<Dipole>::iterator,list<Dipole>::iterator>) const = 0;

  /**
   * For the given performed splitting, dipole chain
   * and dipole taking a recoil, set the next
   * scale.
   */
  virtual void setEvolutionScale(Energy scale,
				 const DipoleSplittingInfo&,
				 DipoleChain&,
				 list<Dipole>::iterator) const = 0;

  /**
   * For the given selected splitting return the evolution scale.
   */
  virtual Energy evolutionScale(const DipoleSplittingInfo&,
				const DipoleSplittingKernel&) const = 0;

  /**
   * Return the maximum pt corresponding to the given
   * evolution scale.
   */
  virtual Energy maxPt(Energy scale,
		       const DipoleSplittingInfo&,
		       const DipoleSplittingKernel&) const = 0;

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
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DipoleEvolutionOrdering> initDipoleEvolutionOrdering;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleEvolutionOrdering & operator=(const DipoleEvolutionOrdering &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleEvolutionOrdering. */
template <>
struct BaseClassTrait<Herwig::DipoleEvolutionOrdering,1> {
  /** Typedef of the first base class of DipoleEvolutionOrdering. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleEvolutionOrdering class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleEvolutionOrdering>
  : public ClassTraitsBase<Herwig::DipoleEvolutionOrdering> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleEvolutionOrdering"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleEvolutionOrdering is implemented. It may also include several, space-separated,
   * libraries if the class DipoleEvolutionOrdering depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleEvolutionOrdering_H */
