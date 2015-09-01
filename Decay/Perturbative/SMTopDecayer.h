// -*- C++ -*-
//
// SMTopDecayer.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SMTopDecayer_H
#define HERWIG_SMTopDecayer_H
//
// This is the declaration of the SMTopDecayer class.
//

#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig/Decay/DecayPhaseSpaceMode.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/Couplings/ShowerAlpha.fh"

namespace Herwig {
  using namespace ThePEG;
  using namespace ThePEG::Helicity;
  
/**
 * \ingroup Decay
 *
 * The SMTopDecayer performs decays of the top quark into
 * the bottom quark and qqbar pairs or to the bottom quark and lepton 
 * neutrino pairs via W boson exchange.
 */
class SMTopDecayer: public DecayIntegrator {

public:

  /**
   * The default constructor.
   */
  SMTopDecayer();

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
  virtual void initializeMECorrection(ShowerTreePtr , double & ,
				      double & );

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

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
   */
  virtual int modeNumber(bool & , tcPDPtr , const tPDVector & ) const {return -1;}

  /**
   * Check if this decayer can perfom the decay for a particular mode.
   * Uses the modeNumber member but can be overridden
   * @param parent The decaying particle
   * @param children The decay products
   */
  virtual bool accept(tcPDPtr parent, const tPDVector & children) const;

  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. As this is the base class this
   * is not implemented.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const Particle & parent,
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
		     const ParticleVector & decay, MEOption meopt) const;

  /**
   * Method to return an object to calculate the 3 (or higher body) partial width
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;
  
  /**
   * The differential three body decay rate with one integral performed.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s  The invariant mass which still needs to be integrate over.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The differential rate \f$\frac{d\Gamma}{ds}\f$
   */
  virtual InvEnergy threeBodydGammads(const int imode, const Energy2 q2,
				      const Energy2 s, const Energy m1,
				      const Energy m2, const Energy m3) const;

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  The integrand for the integrate partial width
   */
  Energy6 dGammaIntegrand(Energy2 mffb2, Energy2 mbf2, Energy mt, Energy mb, 
			  Energy mf, Energy mfb, Energy mw) const;

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
  vector<Lorentz5Momentum> applyHard(const ParticleVector &p,double,double);

  /**
   *  Get the weight for hard emission
   */
  double getHard(double, double);

  /**
   *  This function is auxiliary to the function \f$x_{a}\f$ (hXAB).
   */
  double xgbr(int);

  /**
   *  This function is auxiliary to the function \f$x_{a}\f$ (hXAB).
   */
  double ktr(double,int);

  /**
   *  This function determines \f$x_{a}\f$ as a function of \f$x_{g}\f$ 
   *  and \f$\kappa\f$ where \f$\kappa\f$ pertains to emissions from the 
   *  b.
   */
  double xab(double,double,int);

  /**
   *  This function determines the point (\f$x_{g}\f$) where the condition that 
   *  \f$x_{a}\f$ be real supersedes that due to the external input 
   *  \f$\tilde{\kappa}\f$ where, again, \f$\kappa\f$ pertains to emissions from the 
   *  b.
   */
  double xgbcut(double);

  /**
   *  This function determines the minimum value of \f$x_{a}\f$ 
   *  for a given \f$\tilde{\kappa}\f$ where \f$\kappa\f$ pertains to
   *  emissions from the c.
   */
  double xaccut(double);

  /**
   *  This function is auxiliary to the function \f$x_{g}\f$ (hXGC).
   */
  double z(double,double,int,int); 

  /**
   *  This function determines \f$x_{g}\f$ as a function of \f$x_{a}\f$ 
   *  and \f$\kappa\f$ where \f$\kappa\f$ pertains to emissions from the 
   *  c. It is multivalued, one selects a branch according to the
   *  second to last integer flag (+/-1). The last integer flag
   *  is used to select whether (1) or not (0) you wish to have the 
   *  function for the special case of the full phase space, in which
   *  case the fifth argument \f$\kappa\f$ is irrelevant.
   */
  double xgc(double,double,int,int); 

  /**
   *  This function, \f$x_{g,c=0}^{-1}\f$, returns \f$x_{a}\f$ as a function 
   *  of \f$x_{g}\f$ for the special case of c=0, for emissions from c 
   *  (the b-quark). The third input is \f$\tilde{\kappa}\f$ which pertains 
   *  to emissions from c.
   */
  double xginvc0(double,double); 

  /**
   *  For a given value of \f$x_{g}\f$ this returns the maximum value of \f$x_{a}\f$  
   *  in the dead region.
   */
  double approxDeadMaxxa(double,double,double); 

  /**
   *  For a given value of \f$x_{g}\f$ this returns the maximum value of \f$x_{a}\f$  
   *  in the dead region.
   */
  double approxDeadMinxa(double,double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are in the allowed region, the kinematically accessible phase 
   *  space.
   */
  bool inTheAllowedRegion(double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the approximate dead region.
   */
  bool inTheApproxDeadRegion(double,double,
                                    double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the dead region.
   */
  bool inTheDeadRegion(double,double,
                              double,double); 

  /**
   *  This function returns values of (\f$x_{g}\f$,\f$x_{a}\f$) distributed 
   *  according to \f$\left(1+a-x_{a}\right)^{-1}x_{g}^{-2}\f$ in the 
   *  approximate dead region.  
   */
  double deadRegionxgxa(double,double); 

  /**
   *  This rotation takes a 5-momentum and returns a rotation matrix 
   *  such that it acts on the input 5-momentum so as to
   *  make it point in the +Z direction. Finally it performs a randomn
   *  rotation about the z-axis.
   */
  LorentzRotation rotateToZ(Lorentz5Momentum);

  /**
   *  Full matrix element with a factor of \f$\frac{\alpha_SC_F}{x_g^2\pi}\f$ removed.
   * @param xw The momentum fraction of the W boson
   * @param xg The momentum fraction of the gluon.
   */
  double me(double xw, double xg);

  /**
   *  Access to the strong coupling
   */
  ShowerAlphaPtr coupling() { return _alpha;}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SMTopDecayer> initSMTopDecayer;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SMTopDecayer & operator=(const SMTopDecayer &);
  
  /**
   *Pointer to the W vertex
   */
  AbstractFFVVertexPtr _wvertex;
  
  /**
   * Max weight for integration
   */
  //@{   
  /**
   * Weight \f$W\to q\bar{q}'\f$
   */
  vector<double> _wquarkwgt;
  
  /**
   * Weight \f$W\to \ell \nu\f$
   */
  vector<double> _wleptonwgt;
  //@}

  /**
   *  Pointer to the \f$W^\pm\f$
   */
  PDPtr _wplus;

  /**
   *  Spin density matrix for the decay
   */
  mutable RhoDMatrix _rho;

  /**
   *  1st spinor for the decay
   */
  mutable vector<SpinorWaveFunction   >   _inHalf;

  /**
   *  2nd spinor for the decay
   */
  mutable vector<SpinorWaveFunction   >   _outHalf;

  /**
   *  1st barred spinor for the decay
   */
  mutable vector<SpinorBarWaveFunction>   _inHalfBar;

  /**
   *  2nd barred spinor for the decay
   */
  mutable vector<SpinorBarWaveFunction>   _outHalfBar;

  /**
   *  The mass of the W boson
   */
  Energy _ma;

  /**
   *  The mass of the bottom quark
   */
  Energy _mc;

  /**
   *  The top mass
   */
  Energy _mt;

  /**
   *  The gluon mass.
   */
  Energy _mg;

  /**
   *  The mass ratio for the W.
   */
  double _a;

  /**
   *  The mass ratio for the bottom.
   */
  double _c;

  /**
   *  The mass ratio for the gluon.
   */
  double _g;

  /**
   *  Two times the energy fraction of a.
   */
  double _ktb;

  /**
   *  Two times the energy fraction of the gluon.
   */
  double _ktc;

  /**
   *  Two times the energy fraction of the gluon.
   */
  double _xg;

  /**
   *  Two times the energy fraction of a.
   */
  double _xa;

  /**
   *  Two times the energy fraction of c.
   */
  double _xc;

  /**
   *  This determines the hard matrix element importance 
   *  sampling in _xg. _xg_sampling=2.0 samples as 1/xg^2.
   */
  double _xg_sampling;

  /**
   *  The enhancement factor for initial-state radiation
   */
  double _initialenhance;

  /**
   *  The enhancement factor for final-state radiation
   */
  double _finalenhance;

  /**
   *  This flag determines whether the T2 region in the decay shower
   *  (JHEP12(2003)_045) is populated by the ME correction (true) or
   *  the shower from the decaying particle.
   */
  bool _useMEforT2;

  /**
   *  Pointer to the coupling
   */
  ShowerAlphaPtr _alpha;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SMTopDecayer. */
template <>
struct BaseClassTrait<Herwig::SMTopDecayer,1> {
  /** Typedef of the first base class of SMTopDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SMTopDecayer class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SMTopDecayer>
  : public ClassTraitsBase<Herwig::SMTopDecayer> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::SMTopDecayer"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SMTopDecayer class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwPerturbativeDecay.so"; }
};

/** @endcond */

}

#endif /* HERWIG_SMTopDecayer_H */
