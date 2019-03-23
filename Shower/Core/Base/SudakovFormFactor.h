// -*- C++ -*-
//
// SudakovFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_SudakovFormFactor_H
#define HERWIG_SudakovFormFactor_H
//
// This is the declaration of the SudakovFormFactor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig/Shower/Core/SplittingFunctions/SplittingFunction.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"
#include "Herwig/Shower/Core/SplittingFunctions/SplittingGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ShowerKinematics.fh"
#include "SudakovFormFactor.fh"

namespace Herwig {

using namespace ThePEG;

/**
 *  A typedef for the BeamParticleData
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

/**  \ingroup Shower
 *
 *  This is the definition of the Sudakov form factor class. In general this
 *  is the base class for the implementation of Sudakov form factors in Herwig.
 *  The methods generateNextTimeBranching(), generateNextDecayBranching() and
 *  generateNextSpaceBranching need to be implemented in classes inheriting from this
 *  one.
 *
 *  In addition a number of methods are implemented to assist with the calculation
 *  of the form factor using the veto algorithm in classes inheriting from this one.
 *  
 *  In general the Sudakov form-factor, for final-state radiation, is given
 *  by
 *  \f[\Delta_{ba}(\tilde{q}_{i+1},\tilde{q}_i)=
 *  \exp\left\{
 *     -\int^{\tilde{q}^2_i}_{\tilde{q}^2_{i+1}}
 *     \frac{{\rm d}\tilde{q}^2}{\tilde{q}^2} 
 *      \int\frac{\alpha_S(z,\tilde{q})}{2\pi}
 *      P_{ba}(z,\tilde{q})\Theta(p_T)
 *      \right\}.
 *  \f]
 *  We can solve this to obtain the next value of the scale \f$\tilde{q}_{i+1}\f$
 *  given the previous value \f$\tilde{q}_i\f$
 *  in the following way. First we obtain a simplified form of the integrand
 *  which is greater than or equal to the true integrand for all values of
 *  \f$\tilde{q}\f$.
 *
 *  In practice it is easiest to obtain this over estimate in pieces. The ShowerAlpha
 *  object contains an over estimate for \f$\alpha_S\f$, the splitting function
 *  contains both an over estimate of the spltting function and its integral
 *  which is needed to compute the over estimate of the \f$\tilde{q}\f$ integrand,
 *  together with an over estimate of the limit of the \f$z\f$ integral.
 *
 *  This gives an overestimate of the integrand
 *  \f[g(\tilde{q}^2) = \frac{c}{\tilde{q}^2}, \f]
 *  where because the over estimates are chosen to be independent of \f$\tilde{q}\f$ the 
 *  parameter 
 *  \f[c = \frac{\alpha_{\rm over}}{2\pi}\int^{z_1}_{z_0}P_{\rm over}(z),\f]
 * is a constant independent of \f$\tilde{q}\f$.
 *
 *  The guesst() member can then be used to generate generate the value of 
 *  \f$\tilde{q}^2\f$ according to this result. This is done by solving the Sudakov
 *  form factor, with the over estimates, is equal to a random number 
 *  \f$r\f$ in the interval \f$[0,1]\f$. This gives
 *  \f[\tilde{q}^2_{i+1}=G^{-1}\left[G(\tilde{q}^2_i)+\ln r\right],\f]
 *  where \f$G(\tilde{q}^2)=c\ln(\tilde{q}^2)\f$ is the infinite integral 
 *  of \f$g(\tilde{q}^2)\f$ and \f$G^{-1}(x)=\exp\left(\frac{x}c\right)\f$
 *  is its inverse.
 *  It this case we therefore obtain
 *  \f[\tilde{q}^2_{i+1}=\tilde{q}^2_ir^{\frac1c}.\f]
 *  The value of \f$z\f$ can then be calculated in a similar way
 *  \f[z = I^{-1}\left[I(z_0)+r\left(I(z_1)-I(z_0)\right)\right],\f]
 *  using the guessz() member,
 *  where \f$I=\int P(z){\rm d}z\f$ and \f$I^{-1}\f$ is its inverse.
 *  
 *  The veto algorithm then uses rejection using the ratio of the 
 *  true value to the overestimated one to obtain the original distribution.
 *  This is accomplished using the 
 *  - alphaSVeto()      member for the \f$\alpha_S\f$ veto
 *  - SplittingFnVeto() member for the veto on the value of the splitting function.
 *  in general there must also be a chech that the emission is in the allowed
 *  phase space but this is left to the inheriting classes as it will depend
 *  on the ordering variable.
 *
 *  The Sudakov form factor for the initial-scale shower is different because
 *  it must include the PDF which guides the backward evolution.
 *  It is given by
 *  \f[\Delta_{ba}(\tilde{q}_{i+1},\tilde{q}_i)=
 *  \exp\left\{
 *     -\int^{\tilde{q}^2_i}_{\tilde{q}^2_{i+1}}
 *     \frac{{\rm d}\tilde{q}^2}{\tilde{q}^2} 
 *      \int\frac{\alpha_S(z,\tilde{q})}{2\pi}
 *      P_{ba}(z,\tilde{q})\frac{x'f_a(\frac{x}z,\tilde{q}^2)}{xf_b(x,\tilde{q^2})}
 *      \right\},
 *  \f]
 *  where \f$x\f$ is the fraction of the beam momentum the parton \f$b\f$ had before
 *  the backward evolution.
 *  This can be solve in the same way as for the final-state branching but the constant
 *  becomes
 *  \f[c = \frac{\alpha_{\rm over}}{2\pi}\int^{z_1}_{z_0}P_{\rm over}(z)PDF_{\rm max},\f]
 *  where 
 * \f[PDF_{\rm max}=\max\frac{x'f_a(\frac{x}z,\tilde{q}^2)}{xf_b(x,\tilde{q^2})},\f]
 *  which can be set using an interface.
 *  In addition the PDFVeto() member then is needed to implement the relevant veto.
 *
 *  @see SplittingGenerator
 *  @see SplittingFunction
 *  @see ShowerAlpha
 *  @see \ref SudakovFormFactorInterfaces "The interfaces"
 *  defined for SudakovFormFactor.
 */
class SudakovFormFactor: public Interfaced {

  /**
   *  The SplittingGenerator is a friend to insert the particles in the 
   *  branchings at initialisation
   */
  friend class SplittingGenerator;

public:

  /**
   * The default constructor.
   */
  SudakovFormFactor() : pdfmax_(35.0), pdffactor_(0),
			cutOffOption_(0), a_(0.3), b_(2.3), c_(0.3*GeV),
			kinCutoffScale_( 2.3*GeV ), vgcut_(0.85*GeV),
			vqcut_(0.85*GeV), pTmin_(1.*GeV), pT2min_(ZERO),
			z_( 0.0 ),phi_(0.0), pT_(){}

  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param enhance The radiation enhancement factor
   * defined.
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,
					      const RhoDMatrix & rho,
					      double enhance, double detuning,
					      Energy2 maxQ2)=0;

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextDecayBranching(const Energy startingScale,
					       const Energy stoppingScale,
					       const Energy minmass,
					       const IdList &ids,
					       const RhoDMatrix & rho,
					       double enhance,
					       double detuning)=0;

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns ZERO.
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * defined.
   * @param beam The beam particle
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
					       const IdList &ids,double x,
					       const RhoDMatrix & rho,
					       double enhance,
					       tcBeamPtr beam,
					       double detuning)=0;
  //@}

  /**
   * Generate the azimuthal angle of the branching for forward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiForward(ShowerParticle & particle,const IdList & ids,
				    ShoKinPtr kinematics,
				    const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for backward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiBackward(ShowerParticle & particle,const IdList & ids,
				     ShoKinPtr kinematics,
				     const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for ISR in decays
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiDecay(ShowerParticle & particle,const IdList & ids,
				  ShoKinPtr kinematics,
				  const RhoDMatrix & rho)=0;

  /**
   *  Methods to provide public access to the private member variables
   */
  //@{
  /** 
   * Return the pointer to the SplittingFunction object.
   */
  tSplittingFnPtr splittingFn() const { return splittingFn_; }

  /**
   * Return the pointer to the ShowerAlpha object.
   */
  tShowerAlphaPtr alpha() const { return alpha_; }

  /**
   *  The type of interaction
   */
  inline ShowerInteraction interactionType() const 
  {return splittingFn_->interactionType();}
  //@}

public:

  /**
   *  Methods to access the kinematic variables for the branching
   */
  //@{
  /**
   *  The energy fraction
   */
  double z() const { return z_; }

  /**
   *  The azimuthal angle
   */
  double phi() const { return phi_; }

  /**
   *  The transverse momentum
   */
  Energy pT() const { return pT_; }
  //@}

  /**
   *  Access the maximum weight for the PDF veto
   */
  double pdfMax() const { return pdfmax_;}

  /**
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  virtual Energy calculateScale(double z, Energy pt, IdList ids,unsigned int iopt)=0;

  /**
   *  Method to create the ShowerKinematics object for a final-state branching
   */
  virtual ShoKinPtr createFinalStateBranching(Energy scale,double z,
					      double phi, Energy pt)=0;

  /**
   *  Method to create the ShowerKinematics object for an initial-state branching
   */
  virtual ShoKinPtr createInitialStateBranching(Energy scale,double z,
						double phi, Energy pt)=0;

  /**
   *  Method to create the ShowerKinematics object for a decay branching
   */
  virtual ShoKinPtr createDecayBranching(Energy scale,double z,
					 double phi, Energy pt)=0;

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
  
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

protected:

  /**
   *  Methods to implement the veto algorithm to generate the scale of 
   *  the next branching
   */
  //@{
  /**
   * Value of the energy fraction for the veto algorithm
   * @param iopt The option for calculating z
   * @param ids The PDG codes of the particles in the splitting
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   */
  double guessz (unsigned int iopt, const IdList &ids) const;

  /**
   *  Value of the scale for the veto algorithm
   * @param t1 The starting valoe of the scale
   * @param iopt The option for calculating t 
   * @param ids The PDG codes of the particles in the splitting
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   * @param enhance The radiation enhancement factor
   * @param identical Whether or not the outgoing particles are identical
   */
  Energy2 guesst (Energy2 t1,unsigned int iopt, const IdList &ids,
		  double enhance, bool identical, double detune) const;

  /**
   * Veto on the PDF for the initial-state shower
   * @param t The scale
   * @param x The fraction of the beam momentum
   * @param parton0 Pointer to the particleData for the 
   *                new parent (this is the particle we evolved back to)
   * @param parton1 Pointer to the particleData for the 
   *                original particle
   * @param beam The BeamParticleData object
   */
  bool PDFVeto(const Energy2 t, const double x,
	       const tcPDPtr parton0, const tcPDPtr parton1,
	       tcBeamPtr beam) const;
  /**
   * The PDF veto ratio
   */
  double PDFVetoRatio(const Energy2 t, const double x,
               const tcPDPtr parton0, const tcPDPtr parton1,
               tcBeamPtr beam,double factor) const;

  /**
   *  The veto on the splitting function.
   * @param t The scale
   * @param ids The PDG codes of the particles in the splitting
   * @param mass Whether or not to use the massive splitting functions 
   * @return true if vetoed
   */
  bool SplittingFnVeto(const Energy2 t, 
		       const IdList &ids, 
		       const bool mass,
		       const RhoDMatrix & rho,
		       double detune) const {
    return UseRandom::rnd()>SplittingFnVetoRatio(t,ids,mass,rho,detune);
  }
  
  /**
   * The Splitting function veto ratio
   */
  
  double SplittingFnVetoRatio(const Energy2 t,
			      const IdList &ids,
			      const bool mass,
			      const RhoDMatrix & rho,
			      double detune) const {
    return splittingFn_->ratioP(z_, t, ids,mass,rho)/detune;
  }

  /**
   *  The veto on the coupling constant
   * @param pt2 The value of ther transverse momentum squared, \f$p_T^2\f$.
   * @return true if vetoed
   */
  bool alphaSVeto(Energy2 pt2) const;

  /**
   * The alpha S veto ratio
   */
   
  double alphaSVetoRatio(Energy2 pt2,double factor) const;


  //@}

  /**
   *  Methods to set the kinematic variables for the branching
   */
  //@{
  /**
   *  The energy fraction
   */
  void z(double in) { z_=in; }

  /**
   *  The azimuthal angle
   */
  void phi(double in) { phi_=in; }

  /**
   *  The transverse momentum
   */
  void pT(Energy in) { pT_=in; }
  //@}

  /**
   *  Set/Get the limits on the energy fraction for the splitting
   */
  //@{
  /**
   * Get the limits
   */
  pair<double,double> zLimits() const { return zlimits_;}

  /**
   * Set the limits
   */
  void zLimits(pair<double,double> in) { zlimits_=in; }
  //@}

  /**
   *  Set the particles in the splittings
   */
  void addSplitting(const IdList &);

  /**
   *  Delete the particles in the splittings
   */
  void removeSplitting(const IdList &);

  /**
   *  Access the potential branchings
   */
  const vector<IdList> & particles() const { return particles_; }

public:

  /**
   * @name Methods for the cut-off
   */
  //@{
  /**
   *  The option being used
   */
  unsigned int cutOffOption() const { return cutOffOption_; }

  /**
   *  The kinematic scale
   */
  Energy kinScale() const {return kinCutoffScale_;}

  /**
   * The virtuality cut-off on the gluon \f$Q_g=\frac{\delta-am_q}{b}\f$
   * @param scale The scale \f$\delta\f$
   * @param mq The quark mass \f$m_q\f$.
   */
  Energy kinematicCutOff(Energy scale, Energy mq) const 
  {return max((scale -a_*mq)/b_,c_);}

  /**
   *  The virtualilty cut-off for gluons
   */
  Energy vgCut() const { return vgcut_; }
  
  /**
   *  The virtuality cut-off for everything else
   */
  Energy vqCut() const { return vqcut_; }

  /**
   *  The minimum \f$p_T\f$ for the branching
   */
  Energy pTmin() const { return pTmin_; }
  
  /**
   *  The square of the minimum \f$p_T\f$
   */
  Energy2 pT2min() const { return pT2min_; }

  /**
   *  Calculate the virtual masses for a branchings
   */
  const vector<Energy> & virtualMasses(const IdList & ids);
  //@}

  /**
   *   Set the PDF
   */
  void setPDF(tcPDFPtr pdf, Energy scale) {
    pdf_ = pdf;
    freeze_ = scale;
  }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovFormFactor & operator=(const SudakovFormFactor &) = delete;

private:

  /**
   *  Pointer to the splitting function for this Sudakov form factor
   */
  SplittingFnPtr splittingFn_;

  /**
   *  Pointer to the coupling for this Sudakov form factor
   */
  ShowerAlphaPtr alpha_;

  /**
   * Maximum value of the PDF weight
   */
  double pdfmax_;

  /**
   * List of the particles this Sudakov is used for to aid in setting up
   * interpolation tables if needed
   */
  vector<IdList> particles_;

  /**
   *  Option for the inclusion of a factor \f$1/(1-z)\f$ in the PDF estimate
   */
  unsigned pdffactor_;

private:

  /**
   *  Option for the type of cut-off to be applied
   */
  unsigned int cutOffOption_;

  /**
   *  Parameters for the default Herwig cut-off option, i.e. the parameters for
   *  the \f$Q_g=\max(\frac{\delta-am_q}{b},c)\f$ kinematic cut-off
   */
  //@{
  /**
   *  The \f$a\f$ parameter
   */
  double a_;

  /**
   *  The \f$b\f$ parameter
   */
  double b_;

  /**
   *  The \f$c\f$ parameter
   */
  Energy c_;

  /**
   * Kinematic cutoff used in the parton shower phase space. 
   */
  Energy kinCutoffScale_;
  //@} 

   /**
    *  Parameters for the FORTRAN-like cut-off
    */ 
  //@{
  /**
   *  The virtualilty cut-off for gluons
   */
  Energy vgcut_;
  
  /**
   *  The virtuality cut-off for everything else
   */
  Energy vqcut_;
  //@}

  /**
   *  Parameters for the \f$p_T\f$ cut-off 
   */
  //@{
  /**
   *  The minimum \f$p_T\f$ for the branching
   */
  Energy pTmin_;
  
  /**
   *  The square of the minimum \f$p_T\f$
   */
  Energy2 pT2min_;
  //@}

private:

  /**
   * Member variables to keep the shower kinematics information
   * generated by a call to generateNextTimeBranching or generateNextSpaceBranching
   */
  //@{
  /**
   *  The energy fraction
   */
  double z_;

  /**
   *  The azimuthal angle
   */
  double phi_;

  /**
   *  The transverse momentum
   */
  Energy pT_;
  //@}

  /**
   *  The limits of \f$z\f$ in the splitting
   */
  pair<double,double> zlimits_;

  /**
   *  Stuff for the PDFs
   */
  //@{
  /**
   *  PDf
   */
  tcPDFPtr pdf_;

  /**
   *  Freezing scale
   */
  Energy freeze_;
  //@}


};

}

#endif /* HERWIG_SudakovFormFactor_H */
