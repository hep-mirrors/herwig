// -*- C++ -*-
//
// DipoleChainOrdering.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DipoleChainOrdering_H
#define HERWIG_DipoleChainOrdering_H
//
// This is the declaration of the DipoleChainOrdering class.
//

#include "DipoleEvolutionOrdering.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief DipoleChainOrdering performs ordering on
 * complete colour singlet dipole chains.
 *
 */
class DipoleChainOrdering: public DipoleEvolutionOrdering {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleChainOrdering();

  /**
   * The destructor.
   */
  virtual ~DipoleChainOrdering();
  //@}

public:

  /**
   * For the given dipole and splitting kernel return
   * the hard scale.
   */
  virtual Energy hardScale(tPPtr emitter, tPPtr spectator,
			   double emitterX, double spectatorX,
			   const DipoleSplittingKernel&,
			   const DipoleIndex&) const;

  /**
   * For the given performed splitting, dipole chain
   * and dipoles originating from the splitting, set the next
   * scale.
   */
  virtual void setEvolutionScale(Energy scale,
				 const DipoleSplittingInfo&,
				 DipoleChain&,
				 pair<list<Dipole>::iterator,list<Dipole>::iterator>) const;

  /**
   * For the given performed splitting, dipole chain
   * and dipole taking a recoil, set the next
   * scale.
   */
  virtual void setEvolutionScale(Energy scale,
				 const DipoleSplittingInfo&,
				 DipoleChain&,
				 list<Dipole>::iterator) const;

  /**
   * For the given selected splitting return
   * the evolution scale.
   */
  virtual Energy evolutionScale(const DipoleSplittingInfo& split,
				const DipoleSplittingKernel&) const;

  /**
   * Return the maximum pt corresponding to the given
   * evolution scale.
   */
  virtual Energy maxPt(Energy scale,
		       const DipoleSplittingInfo&,
		       const DipoleSplittingKernel&) const;

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
   * True, if virtuality instead of pt ordering
   * should be performed.
   */
  bool virtualityOrdering;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DipoleChainOrdering> initDipoleChainOrdering;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleChainOrdering & operator=(const DipoleChainOrdering &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DipoleChainOrdering. */
template <>
struct BaseClassTrait<Herwig::DipoleChainOrdering,1> {
  /** Typedef of the first base class of DipoleChainOrdering. */
  typedef Herwig::DipoleEvolutionOrdering NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DipoleChainOrdering class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DipoleChainOrdering>
  : public ClassTraitsBase<Herwig::DipoleChainOrdering> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DipoleChainOrdering"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DipoleChainOrdering is implemented. It may also include several, space-separated,
   * libraries if the class DipoleChainOrdering depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DipoleChainOrdering_H */
