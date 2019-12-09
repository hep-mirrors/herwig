// -*- C++ -*-
//
// HardScaleProfile.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_HardScaleProfile_H
#define Herwig_HardScaleProfile_H
//
// This is the declaration of the HardScaleProfile class.
//

#include "ThePEG/Interface/Interfaced.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief HardScaleProfile is the base class for profile scales. A few
 * standard choices are provided by this implementation.
 *
 */
class HardScaleProfile: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HardScaleProfile();

  /**
   * The destructor.
   */
  virtual ~HardScaleProfile();
  //@}

public:

  /**
   * Return a scale profile towards the hard scale
   */
  virtual double hardScaleProfile(Energy hard, Energy soft) const;

  /**
   * Return true, if this hard scale profile requires an unrestricted
   * radiation phase space.
   */
  virtual bool unrestrictedPhasespace() const;

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

protected:

  /**
   * Enumerate the possible profiles
   */
  enum ProfileTypes {
    theta = 0,
    /** hard theta cut */
    resummation = 1,
    /** `resummation' profile with quadratic interpolation */
    hfact = 2
    /** hfact profile */
  };

  /**
   * A fixed hard scale to be used instead of the hard scale decided
   * for the process in question.
   */
  Energy theFixedHardScale;

  /**
   * A dimensionless width parameter setting the smearing size
   * relative to the hard scale; this may have different
   * interpretations depending on the profile type chosen.
   */
  double theProfileRho;

  /**
   * The profile type to be used
   */
  int theProfileType;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardScaleProfile & operator=(const HardScaleProfile &) = delete;

};

}

#endif /* Herwig_HardScaleProfile_H */
