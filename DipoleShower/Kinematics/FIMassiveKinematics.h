// -*- C++ -*-
//
// FILightKinematics.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_FILightKinematics_H
#define HERWIG_FILightKinematics_H
//
// This is the declaration of the FILightKinematics class.
//

#include "DipoleSplittingKinematics.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup DipoleShower
 * \author Simon Platzer, Martin Stoll
 *
 * \brief FIMassiveKinematics implements massless splittings
 * off a final-initial dipole.
 *
 */
class FIMassiveKinematics: public DipoleSplittingKinematics {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FIMassiveKinematics();

  /**
   * The destructor.
   */
  virtual ~FIMassiveKinematics();
  //@}

public:

  /**
   * Return the boundaries in between the evolution
   * variable random number is to be sampled; the lower
   * cuoff is assumed to correspond to the infrared cutoff.
   */
  virtual pair<double,double> kappaSupport(const DipoleSplittingInfo& dIndex) const;

  /**
   * Return the boundaries in between the momentum
   * fraction random number is to be sampled.
   */
  virtual pair<double,double> xiSupport(const DipoleSplittingInfo& dIndex) const;

  /**
   * Return the dipole scale associated to the
   * given pair of emitter and spectator. This
   * should be the invariant mass or absolute value
   * final/final or initial/initial and the absolute
   * value of the momentum transfer for intial/final or
   * final/initial dipoles.
   */
  virtual Energy dipoleScale(const Lorentz5Momentum& pEmitter,
			     const Lorentz5Momentum& pSpectator) const;

  /**
   * Return the maximum pt for the given dipole scale.
   */
  virtual Energy ptMax(Energy dScale, 
		       double emX, double specX,
		       const DipoleIndex& dIndex,
		       const DipoleSplittingKernel&) const;

  /**
   * Return the maximum virtuality for the given dipole scale.
   */
  virtual Energy QMax(Energy dScale, 
		      double emX, double specX,
		      const DipoleIndex& dIndex) const;

  /**
   * Return the pt given a virtuality.
   */
  virtual Energy PtFromQ(Energy scale, const DipoleSplittingInfo&) const;

  /**
   * Return the virtuality given a pt.
   */
  virtual Energy QFromPt(Energy scale, const DipoleSplittingInfo&) const;

  /**
   * Return the random number associated to
   * the given pt.
   */
  virtual double ptToRandom(Energy pt, Energy dScale,
			    const DipoleIndex& dIndex) const;

  /**
   * Generate splitting variables given three random numbers
   * and the momentum fractions of the emitter and spectator.
   * Return true on success.
   */
  virtual bool generateSplitting(double kappa, double xi, double phi,
				 DipoleSplittingInfo& dIndex);

  /**
   * For the splitting products present in the given dipole splitting
   * info object calculate the kinematics parameters and return the
   * propagator factor.
   */
  virtual InvEnergy2 setKinematics(DipoleSplittingInfo&) const;

  /**
   * For the splitting parameters given in the dipole splitting info
   * object, calculate the phasespace Jacobian times the propagator
   * factor.
   */
  virtual double jacobianTimesPropagator(const DipoleSplittingInfo&,
					 Energy) const;

  /**
   * Generate the full kinematics given emitter and
   * spectator momentum and a previously completeted
   * DipoleSplittingInfo object.
   */
  virtual void generateKinematics(const Lorentz5Momentum& pEmitter,
				  const Lorentz5Momentum& pSpectator,
				  const DipoleSplittingInfo& dInfo);

public:

  /**
   * Triangular / Kallen function
   */
  template <class T>
  inline double rootOfKallen (T a, T b, T c) const {
    return sqrt( a*a + b*b + c*c - 2.*( a*b+a*c+b*c ) ); }

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
  static ClassDescription<FIMassiveKinematics> initFIMassiveKinematics;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FIMassiveKinematics & operator=(const FIMassiveKinematics &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FIMassiveKinematics. */
template <>
struct BaseClassTrait<Herwig::FIMassiveKinematics,1> {
  /** Typedef of the first base class of FIMassiveKinematics. */
  typedef Herwig::DipoleSplittingKinematics NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FIMassiveKinematics class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::FIMassiveKinematics>
  : public ClassTraitsBase<Herwig::FIMassiveKinematics> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::FIMassiveKinematics"; }
  /**
   * The name of a file containing the dynamic library where the class
   * FIMassiveKinematics is implemented. It may also include several, space-separated,
   * libraries if the class FIMassiveKinematics depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwDipoleShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_FIMassiveKinematics_H */
