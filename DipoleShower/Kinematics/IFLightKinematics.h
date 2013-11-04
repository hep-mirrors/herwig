// -*- C++ -*-
//
// IFLightKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_IFLightKinematics_H
#define HERWIG_IFLightKinematics_H
//
// This is the declaration of the IFLightKinematics class.
//

#include "DipoleSplittingKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer
 *
 * \brief IFLightKinematics implements massless splittings
 * off an initial-final dipole.
 *
 * @see \ref IFLightKinematicsInterfaces "The interfaces"
 * defined for IFLightKinematics.
 */
class IFLightKinematics: public DipoleSplittingKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  IFLightKinematics();

  /**
   * The destructor.
   */
  virtual ~IFLightKinematics();
  //@}

public:

  /**
   * Return the maximum pt for the given dipole scale.
   */
  virtual Energy ptMax(Energy dScale, 
		       double emX, double specX,
		       const DipoleIndex&,
		       const DipoleSplittingKernel&) const;

  /**
   * Return the maximum virtuality for the given dipole scale.
   */
  virtual Energy QMax(Energy dScale, 
		      double emX, double specX,
		      const DipoleIndex& dIndex,
		      const DipoleSplittingKernel&) const;

  /**
   * Return the pt given a virtuality.
   */
  virtual Energy PtFromQ(Energy scale, const DipoleSplittingInfo&) const;

  /**
   * Return the virtuality given a pt.
   */
  virtual Energy QFromPt(Energy scale, const DipoleSplittingInfo&) const;

  /**
   * Return the boundaries on the momentum fraction
   */
  virtual pair<double,double> zBoundaries(Energy pt,
					  const DipoleSplittingInfo& dInfo,
					  const DipoleSplittingKernel& split) const;

  /**
   * Generate splitting variables given three random numbers
   * and the momentum fractions of the emitter and spectator.
   * Return true on success.
   */
  virtual bool generateSplitting(double kappa, double xi, double phi,
				 DipoleSplittingInfo& dIndex,
				 const DipoleSplittingKernel&);

  /**
   * Generate the full kinematics given emitter and
   * spectator momentum and a previously completeted
   * DipoleSplittingInfo object.
   */
  virtual void generateKinematics(const Lorentz5Momentum& pEmitter,
				  const Lorentz5Momentum& pSpectator,
				  const DipoleSplittingInfo& dInfo);

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
  static ClassDescription<IFLightKinematics> initIFLightKinematics;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  IFLightKinematics & operator=(const IFLightKinematics &);

private:

  /**
   * Wether or not to choose the `collinear' scheme
   */
  bool theCollinearScheme;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of IFLightKinematics. */
template <>
struct BaseClassTrait<Herwig::IFLightKinematics,1> {
  /** Typedef of the first base class of IFLightKinematics. */
  typedef Herwig::DipoleSplittingKinematics NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the IFLightKinematics class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::IFLightKinematics>
  : public ClassTraitsBase<Herwig::IFLightKinematics> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::IFLightKinematics"; }
  /**
   * The name of a file containing the dynamic library where the class
   * IFLightKinematics is implemented. It may also include several, space-separated,
   * libraries if the class IFLightKinematics depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_IFLightKinematics_H */
