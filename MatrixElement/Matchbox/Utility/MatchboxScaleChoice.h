// -*- C++ -*-
//
// MatchboxScaleChoice.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2012 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef Herwig_MatchboxScaleChoice_H
#define Herwig_MatchboxScaleChoice_H
//
// This is the declaration of the MatchboxScaleChoice class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Handlers/LastXCombInfo.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief MatchboxScaleChoice is the base class for scale choices
 * within Matchbox.
 *
 */
class MatchboxScaleChoice: public HandlerBase, public LastXCombInfo<StandardXComb> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MatchboxScaleChoice();

  /**
   * The destructor.
   */
  virtual ~MatchboxScaleChoice();
  //@}

public:

  /**
   * Clone this scale choice.
   */
  Ptr<MatchboxScaleChoice>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<MatchboxScaleChoice>::ptr>(clone());
  }

  /**
   * Set the XComb object.
   */
  virtual void setXComb(tStdXCombPtr xc) { 
    theLastXComb = xc;
  }

  /**
   * Return the renormalization scale. This default version returns
   * shat.
   */
  virtual Energy2 renormalizationScale() const { 
    return theFixedScale == ZERO ? lastSHat() : sqr(theFixedScale); 
  }

  /**
   * Return the factorization scale. This default version returns
   * shat.
   */
  virtual Energy2 factorizationScale() const { 
    return theFixedScale == ZERO ? lastSHat() : sqr(theFixedScale);
  }

  /**
   * Return the QED renormalization scale. This default version returns
   * the Z mass squared.
   */
  virtual Energy2 renormalizationScaleQED() const { 
    if ( theFixedQEDScale != ZERO )
      return sqr(theFixedQEDScale);
    Energy mZ = getParticleData(ParticleID::Z0)->hardProcessMass();
    return mZ*mZ; 
  }

  /**
   * Return the shower hard scale. This default implementation returns the
   * factorization scale.
   */
  virtual Energy2 showerScale() const {
    return factorizationScale();
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
   * A fixed scale choice. If zero, shat will be used.
   */
  Energy theFixedScale;

  /**
   * A fixed QED scale choice. If zero, shat will be used.
   */
  Energy theFixedQEDScale;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MatchboxScaleChoice & operator=(const MatchboxScaleChoice &);

};

}

#endif /* Herwig_MatchboxScaleChoice_H */
