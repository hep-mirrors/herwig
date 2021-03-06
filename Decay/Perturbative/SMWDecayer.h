// -*- C++ -*-
//
// SMWDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMWDecayer_H
#define HERWIG_SMWDecayer_H
//
// This is the declaration of the SMWDecayer class.
//
#include "Herwig/Decay/PerturbativeDecayer.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig/Decay/PhaseSpaceMode.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Decay
 *
 *  The <code>SMWDecayer</code> is designed to perform the decay of the 
 *  W boson to the Standard Model fermions, including the first order
 *  electroweak corrections.
 *
 * @see PerturbativeDecayer
 * 
 */
class SMWDecayer: public PerturbativeDecayer {

public:

  /**
   * Default constructor.
   */
  SMWDecayer();

public:

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return true;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(RealEmissionProcessPtr , double & ,
				      double & );

  /**
   * Apply the soft matrix element correction
   * @param parent The initial particle in the current branching
   * @param progenitor The progenitor particle of the jet
   * @param fs Whether the emission is initial or final-state
   * @param highestpT The highest pT so far in the shower
   * @param ids ids of the particles produced in the branching
   * @param z The momentum fraction of the branching
   * @param scale the evolution scale of the branching
   * @param pT The transverse momentum of the branching
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(PPtr parent,
				     PPtr progenitor,
				     const bool & fs,
				     const Energy & highestpT,
				     const vector<tcPDPtr> & ids,
				     const double & z,
				     const Energy & scale,
				     const Energy & pT);
  
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

public:

  /**
   * Which of the possible decays is required
   * @param cc Is this mode the charge conjugate
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual int modeNumber(bool & cc, tcPDPtr parent, 
			 const tPDVector & children) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,const tPDVector & children) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param outgoing The particles produced in the decay
   * @param momenta  The momenta of the particles produced in the decay
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(const int ichan,const Particle & part,
	     const tPDVector & outgoing,
	     const vector<Lorentz5Momentum> & momenta,
	     MEOption meopt) const;

  /**
   *   Construct the SpinInfos for the particles produced in the decay
   */
  virtual void constructSpinInfo(const Particle & part,
				 ParticleVector outgoing) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

protected:
  
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}
  
protected:
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   *  Set the \f$\rho\f$ parameter
   */
  void setRho(double);

  /**
   *  Set the \f$\tilde{\kappa}\f$ parameters symmetrically 
   */
  void setKtildeSymm();

  /**
   * Set second \f$\tilde{\kappa}\f$, given the first.
   */
  void setKtilde2();

  /**
   *  Translate the variables from \f$x_q,x_{\bar{q}}\f$ to \f$\tilde{\kappa},z\f$
   */
  //@{
  /**
   *  Calculate \f$z\f$.
   */
  double getZfromX(double, double);

  /**
   *  Calculate \f$\tilde{\kappa}\f$.
   */
  double getKfromX(double, double);
  //@}

  /**
   * Calculate \f$x_{q},x_{\bar{q}}\f$ from \f$\tilde{\kappa},z\f$.
   * @param kt \f$\tilde{\kappa}\f$
   * @param z \f$z\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  void getXXbar(double kt, double z, double & x, double & xbar);

  /**
   *  Soft weight
   */
  //@{
  /**
   *  Soft quark weight calculated from \f$x_{q},x_{\bar{q}}\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  double qWeight(double x, double xbar); 

  /**
   *  Soft antiquark weight calculated from \f$x_{q},x_{\bar{q}}\f$
   * @param x \f$x_{q}\f$
   * @param xbar \f$x_{\bar{q}}\f$
   */
  double qbarWeight(double x, double xbar);

  /**
   * Soft quark weight calculated from \f$\tilde{q},z\f$
   * @param qtilde  \f$\tilde{q}\f$
   * @param z \f$z\f$
   */
  double qWeightX(Energy qtilde, double z);

  /**
   * Soft antiquark weight calculated from \f$\tilde{q},z\f$
   * @param qtilde  \f$\tilde{q}\f$
   * @param z \f$z\f$
   */
  double qbarWeightX(Energy qtilde, double z);
  //@}
  
  /**
   * ????
   */
  double u(double);

  /**
   *  Vector and axial vector parts of the matrix element
   */
  //@{
  /**
   *  Vector part of the matrix element
   */
  double MEV(double, double);

  /**
   *  Axial vector part of the matrix element
   */
  double MEA(double, double);

  /**
   * The matrix element, given \f$x_1\f$, \f$x_2\f$.
   * @param x1 \f$x_1\f$
   * @param x2 \f$x_2\f$
   */
  double PS(double x1, double x2);
  //@}

protected:

  /**
   *  Real emission term, for use in generating the hardest emission
   */
  double calculateRealEmission(double x1, double x2,
			       tcPDVector outgoing,
			       vector<Lorentz5Momentum> momenta,
			       double phi, double muj, double muk,
			       int iemit, bool subtract) const;

  /**
   *  Calculate the ratio between NLO & LO ME
   */
  double meRatio(tcPDVector partons, 
		 vector<Lorentz5Momentum> momenta,
		 unsigned int iemitter,bool subtract) const;

  /**
   *  Calculate matrix element ratio R/B
   */
  virtual double matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				    const ParticleVector & decay3, MEOption meopt,
				    ShowerInteraction inter);
  
  /**
   *  Calculate the LO ME
   */
  double loME(const vector<tcPDPtr> & partons, 
	      const vector<Lorentz5Momentum> & momenta) const;

  /**
   *  Calculate the NLO real emission piece of ME
   */
  InvEnergy2 realME(const vector<tcPDPtr> & partons, 
		    const vector<Lorentz5Momentum> & momenta,
		    ShowerInteraction inter) const;

private:

  /**
   * Private and non-existent assignment operator.
   */
  SMWDecayer & operator=(const SMWDecayer &) = delete;

 private:


  /**
   *  Pointer to the fermion-antifermion W vertex
   */
  AbstractFFVVertexPtr FFWVertex_;

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractFFVVertexPtr FFGVertex_;

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractFFVVertexPtr FFPVertex_;

  /**
   *  Pointer to the fermion-antifermion G vertex
   */
  AbstractVVVVertexPtr WWWVertex_;

  /**
   * maximum weights for the different integrations
   */
  //@{
  /**
   *  Weights for the W to quarks decays.
   */
  vector<double> quarkWeight_;

  /**
   *  Weights for the W to leptons decays.
   */
  vector<double> leptonWeight_;
  //@}

  /**
   *  Spin density matrix for the decay
   */
  mutable RhoDMatrix rho_;

  /**
   *  Polarization vectors for the decay
   */
  mutable vector<VectorWaveFunction> vectors_;

  /**
   *  Spinors for the decay
   */
  mutable vector<SpinorWaveFunction> wave_;

  /**
   *  Barred spinors for the decay
   */
  mutable vector<SpinorBarWaveFunction> wavebar_;

private:

  /**
   * CM energy 
   */
  Energy d_Q_;

  /**
   *  Quark mass
   */
  Energy d_m_;

  /**
   * The rho parameter 
   */
  double d_rho_;

  /**
   * The v parameter
   */
  double d_v_;

  /**
   * The initial kappa-tilde values for radiation from the quark
   */
  double d_kt1_;

  /**
   * The initial kappa-tilde values for radiation from the antiquark
   */
  double d_kt2_;

  /**
   *  Cut-off parameter
   */
  static const double EPS_;

private:

  /**
   *  The colour factor 
   */
  double CF_;

  /**
   *  The W mass
   */
  mutable Energy mW_;


  /**
   *  The strong coupling
   */
  mutable double aS_;

  /**
   * The scale
   */
  mutable Energy2 scale_;

  /**
   *  Stuff for the POWHEG correction
   */
  //@{
  /**
   *  ParticleData object for the gluon
   */
  tcPDPtr gluon_;

  /**
   *  The ParticleData objects for the fermions
   */
  vector<tcPDPtr> partons_;

  /**
   * The fermion momenta
   */
  vector<Lorentz5Momentum> quark_;

  /**
   *  The momentum of the radiated gauge boson
   */
  Lorentz5Momentum gauge_;

  /**
   *  The W boson
   */
  PPtr wboson_;

  /**
   *  W mass squared
   */
  Energy2 mw2_;
  //@}

  /**
   *  Whether ro return the LO or NLO result
   */
  bool NLO_;
};

}


#endif /* HERWIG_SMWDecayer_H */
