// -*- C++ -*-
#ifndef Herwig_PerturbativeDecayer_H
#define Herwig_PerturbativeDecayer_H
//
// This is the declaration of the PerturbativeDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "Herwig/Shower/ShowerAlpha.h"
#include "Herwig/Shower/ShowerInteraction.h"

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
  enum dipoleType {FFa, FFc, IFa, IFc, IFba, IFbc, FFg};

  /**
   *   Phase-space region for an emission (assumes \f$a\to b,c\f$
   */
  enum phaseSpaceRegion {emissionFromB,emissionFromC,emissionFromA1,emissionFromA2,deadZone};

  /**
   *  Type of dipole
   */
  struct DipoleType {

    DipoleType() {}

    DipoleType(dipoleType a, ShowerInteraction b)
      : type(a), interaction(b)
    {}
    
    dipoleType type;

    ShowerInteraction interaction;
  };

public:

  /**
   * The default constructor.
   */
  PerturbativeDecayer() : inter_(ShowerInteraction::QCD),
			  pTmin_(GeV), useMEforT2_(true),
			  C_(5.), ymax_(10.), phaseOpt_(1),
			  pT_(ZERO),mb_(ZERO), e_(0.),
			  s_(0.), e2_(0.), s2_(0.), enhance_(1.)
  {}

  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return No;}

  /**
   *  Member to generate the hardest emission in the POWHEG scheme
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual RealEmissionProcessPtr applyHardMatrixElementCorrection(RealEmissionProcessPtr);

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
   *  Calculate matrix element ratio \f$\frac{M^2}{\alpha_S}\frac{|\overline{\rm{ME}}_3|}{|\overline{\rm{ME}}_2|}\f$
   */
  virtual double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				    const ParticleVector & decay3, MEOption meopt,
				    ShowerInteraction inter);

  /**
   *  Work out the type of process
   */
  bool identifyDipoles(vector<DipoleType> & dipoles,
		       PPtr & aProgenitor,
		       PPtr & bProgenitor,
		       PPtr & cProgenitor,
		       ShowerInteraction inter) const;

  /**
   *  Coupling for the generation of hard QCD radiation
   */
  ShowerAlphaPtr alphaS() {return alphaS_;}

  /**
   *  Coupling for the generation of hard QED radiation
   */
  ShowerAlphaPtr alphaEM() {return alphaEM_;}

  /**
   *  Return the momenta including the hard emission
   */
  vector<Lorentz5Momentum> hardMomenta(PPtr in, PPtr emitter, 
  				       PPtr spectator, 
  				       const vector<DipoleType>  & dipoles,
				       int i, bool inDeadZone);

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
   * Return dipole corresponding to the DipoleType dipoleId
   */
  pair<double,double> calculateDipole(const DipoleType & dipoleId,
				      const Particle & inpart,
				      const ParticleVector & decay3);

  /**
   * Return contribution to dipole that depends on the spin of the emitter
   */
  double dipoleSpinFactor(tcPDPtr emitter, double z);

  /**
   *  Return the colour coefficient of the dipole
   */
  double colourCoeff(tcPDPtr emitter, tcPDPtr spectator,
		     tcPDPtr other, DipoleType dipole);
  
  /**
   * Set up the colour lines
   */
  void getColourLines(RealEmissionProcessPtr real);

  /**
   *  Generate a hard emission
   */
  RealEmissionProcessPtr getHardEvent(RealEmissionProcessPtr born,
				      bool inDeadZone,
				      ShowerInteraction inter);

  /**
   *  Is the \f$x_g,x_s\f$ point in the dead-zone for all the dipoles
   */
  bool inTotalDeadZone(double xg, double xs,
		       const vector<DipoleType>  & dipoles,
		       int i);

  /**
   *  Is the \f$x_g,x_a\f$ point in the dead-zone for an initial-final colour connection
   */
  phaseSpaceRegion inInitialFinalDeadZone(double xg, double xa, double a, double c) const;

  /**
   *  Is the \f$x_b,x_c\f$ point in the dead-zone for a final-final colour connection
   */
  phaseSpaceRegion inFinalFinalDeadZone(double xb, double xc, double b, double c) const;

  /**
   *  For me corrections use the shower or me for the T2 region
   */
  bool useMEforT2() const {return useMEforT2_;}

protected:

  /**
   *  Access to the kinematics for inheriting classes
   */
  //@{
  /**
   *  Transverse momentum of the emission
   */
  const Energy & pT() const { return pT_;}

  /**
   *  Mass of decaying particle
   */
  const Energy & mb() const {return mb_;}

  /**
   *  Reduced mass of emitter child particle
   */
  const double & e() const {return e_;}

  /**
   * Reduced mass of spectator child particle
   */
  const double & s() const {return s_;}

  /**
   *  Reduced mass of emitter child particle squared
   */
  const double & e2() const {return e2_;}

  /**
   * Reduced mass of spectator child particle squared
   */
  const double & s2() const {return s2_;}
  //@}
 
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  PerturbativeDecayer & operator=(const PerturbativeDecayer &) = delete;

private:

  /**
   *  Members for the generation of the hard radiation
   */
  //@{
  /**
   *  Which types of radiation to generate
   */
  ShowerInteraction inter_;
  
  /**
   *  Coupling for the generation of hard QCD radiation
   */
  ShowerAlphaPtr alphaS_;
  
  /**
   *  Coupling for the generation of hard QED radiation
   */
  ShowerAlphaPtr alphaEM_;

  /**
   *  Minimum \f$p_T\f$
   */
  Energy pTmin_;

  /**
   *  This flag determines whether the T2 region in the decay shower
   *  (JHEP12(2003)_045) is populated by the ME correction (true) or
   *  the shower from the decaying particle.
   */
  bool useMEforT2_;

  /**
   *   Prefactor for the sampling
   */
  double C_;

  /**
   *   Maximum value for y
   */
  double ymax_;

  /**
   *  Option for phase-space sampling
   */
  unsigned int phaseOpt_;
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

  /**
   *  Enhancement prefactor for special cases
   */
  mutable double enhance_;
  //@}
};

}

#endif /* Herwig_PerturbativeDecayer_H */
