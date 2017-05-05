// -*- C++ -*-
//
// SMZDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMZDecayer_H
#define HERWIG_SMZDecayer_H
//
// This is the declaration of the SMZDecayer class.
//
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.fh"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/** \ingroup Decay
 *
 *  The <code>SMZDecayer</code> is designed to perform the decay of the 
 *  Z boson to the Standard Model fermions. In principle it can also
 *  be used for these decays in any model.
 *
 * @see DecayIntegrator
 * 
 */
class SMZDecayer: public DecayIntegrator {

public:

  /**
   * Default constructor.
   */
  SMZDecayer();

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
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual RealEmissionProcessPtr applyHardMatrixElementCorrection(RealEmissionProcessPtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br);
  //@}

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
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @param meopt Option for the calculation of the matrix element
   * @return The matrix element squared for the phase-space configuration.
   */
  virtual double me2(const int ichan, const Particle & part,
		     const ParticleVector & decay,MEOption meopt) const;

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   */
  virtual void dataBaseOutput(ofstream & os,bool header) const;

  /**
   *  Members for the generation of QED radiation in the decays
   */
  //@{
  /**
   *  The one-loop virtual correction.
   * @param imode The mode required.
   * @param part  The decaying particle.
   * @param products The decay products including the radiated photon.
   * @return Whether the correction is implemented
   */
  virtual double oneLoopVirtualME(unsigned int imode,
				  const Particle & part, 
				  const ParticleVector & products);
  
  /**
   *  The real emission matrix element
   * @param imode The mode required
   * @param part  The decaying particle
   * @param products The decay products including the radiated photon
   * @param iemitter The particle which emitted the photon
   * @param ctheta   The cosine of the polar angle between the photon and the
   *                 emitter
   * @param stheta   The sine of the polar angle between the photon and the
   *                 emitter 
   * @param rot1 Rotation from rest frame to frame for real emission
   * @param rot2 Rotation to place emitting particle along z
   */
  virtual InvEnergy2 realEmissionME(unsigned int imode,
				    const Particle & part, 
				    ParticleVector & products,
				    unsigned int iemitter,
				    double ctheta, double stheta,
				    const LorentzRotation & rot1,
				    const LorentzRotation & rot2);
  //@}

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
   *  Apply the hard matrix element
   */
  vector<Lorentz5Momentum> applyHard(const ParticleVector &p);

  /**
   *  Get the weight for hard emission
   */
  double getHard(double &, double &);

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

  /**
   *  Access to the strong coupling
   */
  ShowerAlphaPtr alphaS() const {return alpha_;}
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SMZDecayer> initSMZDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  SMZDecayer & operator=(const SMZDecayer &);

private:

  /**
   * Pointer to the Z vertex
   */
  FFVVertexPtr FFZvertex_;

  /**
   * Pointer to the photon vertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   * maximum weights for the different integrations
   */
  //@{
  /**
   *  Weights for the Z to quarks decays.
   */
  vector<double> quarkWeight_;

  /**
   *  Weights for the Z to leptons decays.
   */
  vector<double> leptonWeight_;
  //@}

  /**
   *  Spin density matrix for the decay
   */
  mutable RhoDMatrix _rho;

  /**
   *  Polarization vectors for the decay
   */
  mutable vector<VectorWaveFunction> _vectors;

  /**
   *  Spinors for the decay
   */
  mutable vector<SpinorWaveFunction> _wave;

  /**
   *  Barred spinors for the decay
   */
  mutable vector<SpinorBarWaveFunction> _wavebar;

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

  /**
   *  Pointer to the coupling
   */
  ShowerAlphaPtr alpha_;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SMZDecayer.
 */
template <>
 struct BaseClassTrait<Herwig::SMZDecayer,1> {
    /** Typedef of the base class of SMZDecayer. */
   typedef Herwig::DecayIntegrator NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
 struct ClassTraits<Herwig::SMZDecayer>
  : public ClassTraitsBase<Herwig::SMZDecayer> {
   /** Return the class name.*/
  static string className() { return "Herwig::SMZDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwPerturbativeDecay.so"; }

};

/** @endcond */

}

#endif /* HERWIG_SMZDecayer_H */
