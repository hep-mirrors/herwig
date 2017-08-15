// -*- C++ -*-
#ifndef Herwig_PerturbativeDecayer_H
#define Herwig_PerturbativeDecayer_H
//
// This is the declaration of the PerturbativeDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The PerturbativeDecayer class is the base class for perturbative decays in
 * Herwig and implements the functuality for the POWHEG corrections
 *
 * @see \ref PerturbativeDecayerInterfaces "The interfaces"
 * defined for PerturbativeDecayer.
 */
class PerturbativeDecayer: public DecayIntegrator {

protected:

  /**
   * Type of dipole
   */
  enum dipoleType {FFa, FFc, IFa, IFc, IFba, IFbc};

public:

  /**
   * The default constructor.
   */
  PerturbativeDecayer() : pTmin_(GeV), pT_(ZERO), mb_(ZERO),
			  e_(0.), s_(0.), e2_(0.), s2_(0.)
  {}

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}

  /**
   *  Member to generate the hardest emission in the POWHEG scheme
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);

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

  /**
   *  Three-body matrix element including additional QCD radiation
   */
  virtual double threeBodyME(const int , const Particle & inpart,
			     const ParticleVector & decay, MEOption meopt);

  /**
   *  Calculate matrix element ratio R/B
   */
  virtual double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				    const ParticleVector & decay3, MEOption meopt);

  /**
   *  Work out the type of process
   */
  bool identifyDipoles(vector<dipoleType> & dipoles,
		       PPtr & aProgenitor,
		       PPtr & bProgenitor,
		       PPtr & cProgenitor) const;

  /**
   *  Coupling for the generation of hard radiation
   */
  ShowerAlphaPtr coupling() {return coupling_;}

  /**
   *  Return the momenta including the hard emission
   */
  vector<Lorentz5Momentum> hardMomenta(PPtr in, PPtr emitter, 
  				       PPtr spectator, 
  				       const vector<dipoleType>  & dipoles, int i);

  /**
   *  Calculate momenta of all the particles
   */
  bool calcMomenta(int j, Energy pT, double y, double phi, double& xg, 
		   double& xs, double& xe, double& xe_z, 
		   vector<Lorentz5Momentum>& particleMomenta);

  /**
   *  Check the calculated momenta are physical
   */
  bool psCheck(const double xg, const double xs);

  /**
   * Return dipole corresponding to the dipoleType dipoleId
   */
  InvEnergy2 calculateDipole(const dipoleType & dipoleId,   const Particle & inpart,
			     const ParticleVector & decay3, const dipoleType & emittingDipole);

  /**
   * Return contribution to dipole that depends on the spin of the emitter
   */
  double dipoleSpinFactor(const PPtr & emitter, double z);

  /**
   *  Return the colour coefficient of the dipole
   */
  double colourCoeff(const PDT::Colour emitter, const PDT::Colour spectator,
		     const PDT::Colour other);
  
  /**
   * Set up the colour lines
   */
  void getColourLines(RealEmissionProcessPtr real);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PerturbativeDecayer & operator=(const PerturbativeDecayer &);

private:

  /**
   *  Members for the generation of the hard radiation
   */
  //@{
  /**
   *  Coupling for the generation of hard radiation
   */
  ShowerAlphaPtr coupling_;

  /**
   *  Minimum \f$p_T\f$
   */
  Energy pTmin_;
  //@}

private:

  /**
   *   Mmeber variables for the kinematics of the hard emission
   */
  //@{
  /**
   *  Transverse momentum of the emission
   */
  Energy pT_;

  /**
   *  Mass of decaying particle
   */
  Energy mb_;

  /**
   *  Reduced mass of emitter child particle
   */
  double e_;

  /**
   * Reduced mass of spectator child particle
   */
  double s_;

  /**
   *  Reduced mass of emitter child particle squared
   */
  double e2_;

  /**
   * Reduced mass of spectator child particle squared
   */
  double s2_;
  //@}
};

}

#endif /* Herwig_PerturbativeDecayer_H */
