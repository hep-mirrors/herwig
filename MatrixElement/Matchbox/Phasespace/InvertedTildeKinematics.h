// -*- C++ -*-
//
// InvertedTildeKinematics.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_InvertedTildeKinematics_H
#define HERWIG_InvertedTildeKinematics_H
//
// This is the declaration of the InvertedTildeKinematics class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/MatrixElement/Matchbox/Dipoles/SubtractionDipole.h"

namespace Herwig {

using namespace ThePEG;

/**
 * \ingroup Matchbox
 * \author Simon Platzer
 *
 * \brief InvertedTildeKinematics is the base class for the inverted 'tilde'
 * kinematics being used for subtraction terms in the
 * formalism of Catani and Seymour.
 *
 */
class InvertedTildeKinematics: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  InvertedTildeKinematics();

  /**
   * The destructor.
   */
  virtual ~InvertedTildeKinematics();
  //@}

public:

  /** @name Access to kinematic quantities. */
  //@{
  /**
   * Return the momentum of the emitter in the real emission process
   */
  const Lorentz5Momentum& realEmitterMomentum() const { return theRealEmitterMomentum; }

  /**
   * Return the momentum of the emission in the real emission process
   */
  const Lorentz5Momentum& realEmissionMomentum() const { return theRealEmissionMomentum; }

  /**
   * Return the momentum of the spectator in the real emission process
   */
  const Lorentz5Momentum& realSpectatorMomentum() const { return theRealSpectatorMomentum; }

  /**
   * Return the momentum of the emitter in the underlying Born process
   */
  const Lorentz5Momentum& bornEmitterMomentum() const { 
    return theBornXComb->meMomenta()[theDipole->bornEmitter()];
  }

  /**
   * Return the momentum of the spectator in the underlying Born process
   */
  const Lorentz5Momentum& bornSpectatorMomentum() const { 
    return theBornXComb->meMomenta()[theDipole->bornSpectator()];
  }

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

  /**
   * Return the vector of dimensionless variables calculated
   */
  const vector<double>& subtractionParameters() const { return theDipole->subtractionParameters(); }

  /**
   * Return true, if this InvertedTildeKinematics object needs to transform
   * all other particles in the process except the emitter, emission and spectator
   */
  virtual bool doesTransform() const { return false; }

  /**
   * If this InvertedTildeKinematics object needs to transform all other particles
   * in the process except the emitter, emission and spectator, return the transformed
   * momentum.
   */
  virtual Lorentz5Momentum transform(const Lorentz5Momentum& p) const { return p; }

  /**
   * Return the centre of mass energy for the underlying Born configuration
   */
  Energy2 sHat() const { return theBornXComb->lastSHat(); }
  //@}

public:

  /**
   * Clone this object
   */
  Ptr<InvertedTildeKinematics>::ptr cloneMe() const {
    return dynamic_ptr_cast<Ptr<InvertedTildeKinematics>::ptr>(clone());
  }

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
   * Return the real xcomb
   */
  tcStdXCombPtr realXComb() const { return theRealXComb; }

  /**
   * Return the Born xcomb
   */
  tcStdXCombPtr bornXComb() const { return theBornXComb; }

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
   * Return the number of random numbers needed to generate
   * a real emission configuration off the underlying Born
   * configuration.
   */
  virtual int nDimRadiation() const { return 3; }

  /**
   * Perform the mapping of the tilde kinematics for the
   * last selected process and store all dimensionless
   * variables in the subtractionParameters() vector.
   * Return false, if the calculation of the real
   * kinematics was impossible for the selected configuration
   * and true on success.
   */
  virtual bool doMap(const double *) = 0;

  /**
   * Set an optional cutoff on the emission's
   * transverse momentum.
   */
  void ptCut(Energy pt) { thePtCut = pt; }
  
  /**
   * Return the optional cutoff on the emission's
   * transverse momentum.
   */
  Energy ptCut() const { return thePtCut; }

  /**
   * Return the random number index
   * corresponding to the evolution variable.
   */
  virtual int evolutionVariable() const { return 0; }

  /**
   * Return the cutoff on the evolution
   * random number corresponding to the pt cut.
   */
  virtual double evolutionCutoff() const { return 0.0; }

  /**
   * Return the pt associated to the last generated splitting.
   */
  virtual Energy lastPt() const = 0;

  /**
   * Return the momentum fraction associated to the last splitting.
   */
  virtual double lastZ() const = 0;

  /**
   * Return the relevant dipole scale
   */
  virtual Energy lastScale() const;

  /**
   * Return the upper bound on pt
   */
  virtual Energy ptMax() const = 0;

  /**
   * Given a pt and a hard pt, return the boundaries on z; if the hard
   * pt is zero, ptMax() will be used.
   */
  virtual pair<double,double> zBounds(Energy pt, Energy hardPt = ZERO) const = 0;

  /**
   * Generate pt and z
   */
  virtual pair<Energy,double> generatePtZ(double& jac, const double * r,
  					  double power=1., vector<double>* values = NULL) const;

  /**
   * Return the single particle phase space weight in units
   * of sHat() for the last selected configuration.
   */
  double jacobian() const { return theJacobian; }

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
   * Access the momentum of the emitter in the real emission process
   */
  Lorentz5Momentum& realEmitterMomentum() { return theRealEmitterMomentum; }

  /**
   * Access the momentum of the emission in the real emission process
   */
  Lorentz5Momentum& realEmissionMomentum() { return theRealEmissionMomentum; }

  /**
   * Access the momentum of the spectator in the real emission process
   */
  Lorentz5Momentum& realSpectatorMomentum() { return theRealSpectatorMomentum; }

  /**
   * Access the vector of dimensionless variables calculated
   */
  vector<double>& subtractionParameters() { return theDipole->subtractionParameters(); }

  /**
   * Set the single particle phase space weight in units
   * of sHat() for the last selected configuration.
   */
  void jacobian(double w) { theJacobian = w; }

  /**
   * Calculate a transverse momentum for the given momenta,
   * invariant pt and azimuth.
   */
  Lorentz5Momentum getKt(const Lorentz5Momentum& p1,
			 const Lorentz5Momentum& p2,
			 Energy pt,
			 double phi,
			 bool spacelike = false) const;

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
   * The last dipole this InvertedTildeKinematics has been selected for
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
   * The momentum of the emitter in the real emission process
   */
  Lorentz5Momentum theRealEmitterMomentum;

  /**
   * The momentum of the emission in the real emission process
   */
  Lorentz5Momentum theRealEmissionMomentum;

  /**
   * The momentum of the spectator in the real emission process
   */
  Lorentz5Momentum theRealSpectatorMomentum;

  /**
   * Return the single particle phase space weight in units
   * of sHat() for the last selected configuration.
   */
  double theJacobian;

  /**
   * The optional cutoff on the emission's
   * transverse momentum.
   */
  Energy thePtCut;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  InvertedTildeKinematics & operator=(const InvertedTildeKinematics &);

};

}

#endif /* HERWIG_InvertedTildeKinematics_H */
