// -*- C++ -*-
//
// TildeKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_TildeKinematics_H
#define HERWIG_TildeKinematics_H
//
// This is the declaration of the TildeKinematics class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief TildeKinematics is the base class for the 'tilde'
 * kinematics being used for subtraction terms in the
 * formalism of Catani and Seymour.
 *
 */
class TildeKinematics: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  TildeKinematics();

  /**
   * The destructor.
   */
  virtual ~TildeKinematics();
  //@}

public:

  /**
   * Clone this object
   */
  Ptr<TildeKinematics>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<TildeKinematics>::ptr>(clone());
  }

  /** @name Access to kinematic quantities. */
  //@{
  /**
   * Return the momentum of the emitter in the real emission process
   */
  const Lorentz5Momentum& realEmitterMomentum() const {
    return theRealXComb->meMomenta()[theDipole->realEmitter()];
  }

  /**
   * Return the momentum of the emission in the real emission process
   */
  const Lorentz5Momentum& realEmissionMomentum() const {
    return theRealXComb->meMomenta()[theDipole->realEmission()];
  }

  /**
   * Return the momentum of the spectator in the real emission process
   */
  const Lorentz5Momentum& realSpectatorMomentum() const {
    return theRealXComb->meMomenta()[theDipole->realSpectator()];
  }

  /**
   * Return the momentum of the emitter in the underlying Born process
   */
  const Lorentz5Momentum& bornEmitterMomentum() const { return theBornEmitterMomentum; }

  /**
   * Return the momentum of the spectator in the underlying Born process
   */
  const Lorentz5Momentum& bornSpectatorMomentum() const { return theBornSpectatorMomentum; }

  /**
   * Return the vector of dimensionless variables calculated
   */
  const vector<double>& subtractionParameters() const { return theDipole->subtractionParameters(); }

  /**
   * Return true, if this TildeKinematics object needs to transform
   * all other particles in the process except the emitter and spectator
   */
  virtual bool doesTransform() const { return false; }

  /**
   * If this TildeKinematics object needs to transform all other particles
   * in the process except the emitter and spectator, return the transformed
   * momentum.
   */
  virtual Lorentz5Momentum transform(const Lorentz5Momentum& p) const { return p; }
  //@}

  /**
   * If this tilde kinematics is implementing a mapping different from
   * the baseline dipole mapping, determine the relevant shower
   * parameters and check for phase space boundaries. Note that real
   * emission kinematics only are available at this stage.
   */
  virtual void getShowerVariables() const {}

  /**
   * If this tilde kinematics is implementing a mapping different from
   * the baseline dipole mapping, return the ratio of phase space
   * factorization Jacobians for this and the nominal dipole
   * mapping. This is used for matching subtractions.
   */
  virtual double jacobianRatio() const { return 1.; }

public:

  /** @name Access to process data. */
  //@{
  /**
   * Prepare given a dipole, and XCombs describing the real emission
   * and underlying Born processes, respectively.
   */
  void prepare(tcStdXCombPtr newRealXComb,
	       tcStdXCombPtr newBornXComb) {
    theRealXComb = newRealXComb; theBornXComb = newBornXComb;
  }

  /**
   * Set the current dipole
   */
  void dipole(Ptr<SubtractionDipole>::tptr dip) { theDipole = dip; }

  /**
   * Return the current dipole
   */
  Ptr<SubtractionDipole>::tptr dipole() { return theDipole; }

  /**
   * Return the current dipole
   */
  Ptr<SubtractionDipole>::tcptr dipole() const { return theDipole; }

  /**
   * Perform the mapping to the tilde kinematics for the
   * last selected process and store all dimensionless
   * variables in the subtractionParameters() vector.
   * Return false, if the calculation of the tilde
   * kinematics was impossible for the selected configuration
   * and true on success.
   */
  virtual bool doMap() = 0;

  /**
   * Return the pt associated to the last merged splitting.
   */
  virtual Energy lastPt() const = 0;
  
  
  /**
   * Return the pt associated to emitter emission and sppectator momentum.
   */
  virtual Energy lastPt(Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum) const =0 ;


  /**
   * Given a pt and a hard pt, return the boundaries on z; 
   */
  virtual pair<double,double> zBounds(Energy pt, Energy hardPt ) const = 0;
  
  
  /**
   * Return the momentum fraction associated to the last splitting.
   */
  virtual double lastZ() const = 0;

  /**
   * Return the relevant dipole scale
   */
  virtual Energy lastScale() const;
  
  virtual bool aboveAlpha() const {
	cerr<<"only implemented for light kinematics";
        assert(false);
	return false;
  }

  /**
   * Return the particle type of the emitter in the real emission process
   */
  cPDPtr realEmitterData() const { 
    return 
      (theDipole && theRealXComb) ? 
      theRealXComb->mePartonData()[theDipole->realEmitter()] :
      cPDPtr();
  }

  /**
   * Return the particle type of the emission in the real emission process
   */
  cPDPtr realEmissionData() const { 
    return 
      (theDipole && theRealXComb) ? 
      theRealXComb->mePartonData()[theDipole->realEmission()] :
      cPDPtr();
  }

  /**
   * Return the particle type of the spectator in the real emission process
   */
  cPDPtr realSpectatorData() const { 
    return 
      (theDipole && theRealXComb) ? 
      theRealXComb->mePartonData()[theDipole->realSpectator()] :
      cPDPtr();
  }

  /**
   * Return the particle type of the emitter in the underlying Born process
   */
  cPDPtr bornEmitterData() const { 
    return 
      (theDipole && theBornXComb) ? 
      theBornXComb->mePartonData()[theDipole->bornEmitter()] :
      cPDPtr();
  }

  /**
   * Return the particle type of the spectator in the underlying Born process
   */
  cPDPtr bornSpectatorData() const { 
    return 
      (theDipole && theBornXComb) ? 
      theBornXComb->mePartonData()[theDipole->bornSpectator()] :
      cPDPtr();
  }
  //@}

protected:

  /**
   * Access the momentum of the emitter in the underlying Born process
   */
  Lorentz5Momentum& bornEmitterMomentum() { return theBornEmitterMomentum; }

  /**
   * Access the momentum of the spectator in the underlying Born process
   */
  Lorentz5Momentum& bornSpectatorMomentum() { return theBornSpectatorMomentum; }

  /**
   * Access the vector of dimensionless variables calculated
   */
  vector<double>& subtractionParameters() { return theDipole->subtractionParameters(); }

public: 
  /**
   * Return the momentum fraction of the emitter
   */
  double emitterX() const {
    return
    theDipole->bornEmitter() == 0 ?
    theBornXComb->lastX1() :
    theBornXComb->lastX2();
  }
  
  
  /**
   * Return the momentum fraction of the spectator
   */
  double spectatorX() const {
    return
    theDipole->bornSpectator() == 0 ?
    theBornXComb->lastX1() :
    theBornXComb->lastX2();
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


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}


private:

  /**
   * The last dipole this TildeKinematics has been selected for
   */
  Ptr<SubtractionDipole>::tptr theDipole;

  /**
   * The XComb object describing the real emission process
   */
  tcStdXCombPtr theRealXComb;

  /**
   * The XComb object describing the underlying Born process
   */
  tcStdXCombPtr theBornXComb;

  /**
   * The momentum of the emitter in the underlying Born process
   */
  Lorentz5Momentum theBornEmitterMomentum;

  /**
   * The momentum of the spectator in the underlying Born process
   */
  Lorentz5Momentum theBornSpectatorMomentum;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TildeKinematics & operator=(const TildeKinematics &) = delete;

};

}

#endif /* HERWIG_TildeKinematics_H */
