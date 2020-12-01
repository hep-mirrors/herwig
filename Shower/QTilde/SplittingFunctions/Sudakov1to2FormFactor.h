// -*- C++ -*-
//
// Sudakov1to2FormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_Sudakov1to2FormFactor_H
#define HERWIG_Sudakov1to2FormFactor_H
//
// This is the declaration of the Sudakov1to2FormFactor class.
//

#include "SudakovFormFactor.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/SplittingGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "Sudakov1to2FormFactor.fh"
#include "SudakovCutOff.h"
#include "Herwig/Decay/DecayMatrixElement.h"


namespace Herwig {

using namespace ThePEG;

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
 *  The guesstz() member can then be used to generate generate the value of 
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
 *  using the guesstz() member,
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
 *
 *  This is an abstract class which defines the common interface
 *  for all \f$1\to2\f$ splitting functions, for both initial-state
 *  and final-state radiation. 
 *
 *  The SplittingFunction class contains a number of purely virtual members
 *  which must be implemented in the inheriting classes. The class also stores
 *  the interaction type of the spltting function.
 *
 *  The inheriting classes need to specific the splitting function 
 *  \f$P(z,2p_j\cdot p_k)\f$, in terms of the energy fraction \f$z\f$ and
 *  the evolution scale. In order to allow the splitting functions to be used
 *  with different choices of evolution functions the scale is given by
 * \f[2p_j\cdot p_k=(p_j+p_k)^2-m_{jk}^2=Q^2-(p_j+p_k)^2=z(1-z)\tilde{q}^2=
 *   \frac{p_T^2}{z(1-z)}-m_{jk}^2+\frac{m_j^2}{z}+\frac{m_k^2}{1-z},\f]
 *  where \f$Q^2\f$ is the virtuality of the branching particle,
 *  $p_T$ is the relative transverse momentum of the branching products and
 *  \f$\tilde{q}^2\f$ is the angular variable described in hep-ph/0310083.
 *
 *  In addition an overestimate of the 
 *  splitting function, \f$P_{\rm over}(z)\f$ which only depends upon \f$z\f$, 
 *  the integral, inverse of the integral for this overestimate and
 *  ratio of the true splitting function to the overestimate must be provided
 *  as they are necessary for the veto alogrithm used to implement the evolution.
 *
 *  In addition the PDFVeto() member then is needed to implement the relevant veto.
 *
 *  @see SplittingGenerator
 *  @see ShowerAlpha
 *  @see \ref Sudakov1to2FormFactorInterfaces "The interfaces"
 *  defined for Sudakov1to2FormFactor.
 */
class Sudakov1to2FormFactor: public SudakovFormFactor {

public:

  /**
   * The default constructor.
   */
  Sudakov1to2FormFactor() : 
			z_( 0.0 ),phi_(0.0), pT_(), scaleChoice_(2),
			strictAO_(true), colourFactor_(-1.)
  {}

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
  ShoKinPtr generateNextTimeBranching(const Energy startingScale,
				      const IdList &ids,
				      const RhoDMatrix & rho,
				      double enhance, double detuning);

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
  ShoKinPtr generateNextDecayBranching(const Energy startingScale,
				       const Energy stoppingScale,
				       const Energy minmass,
				       const IdList &ids,
				       const RhoDMatrix & rho,
				       double enhance,
				       double detuning);

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
  ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
				       const IdList &ids,double x,
				       const RhoDMatrix & rho,
				       double enhance,
				       tcBeamPtr beam,
				       double detuning);
  //@}

  /**
   * Generate the azimuthal angle of the branching for forward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  double generatePhiForward(ShowerParticle & particle,const IdList & ids,
			    ShoKinPtr kinematics,
			    const RhoDMatrix & rho);

  /**
   *  Generate the azimuthal angle of the branching for backward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  double generatePhiBackward(ShowerParticle & particle,const IdList & ids,
			     ShoKinPtr kinematics,
			     const RhoDMatrix & rho);

  /**
   *  Generate the azimuthal angle of the branching for ISR in decays
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  double generatePhiDecay(ShowerParticle & particle,const IdList & ids,
			  ShoKinPtr kinematics,
			  const RhoDMatrix & rho);
  //@}

public:

  /**
   *   Methods to return the splitting function.
   */
  //@{
  /**
   * Purely virtual method which should return the exact value of the splitting function,
   * \f$P\f$ evaluated in terms of the energy fraction, \f$z\f$, and the evolution scale 
   \f$\tilde{q}^2\f$.
   * @param z   The energy fraction.
   * @param t   The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double P(const double z, const Energy2 t, const IdList & ids,
		   const bool mass, const RhoDMatrix & rho) const = 0;

  /**
   * Purely virtual method which should return
   * an overestimate of the splitting function,
   * \f$P_{\rm over}\f$ such that the result \f$P_{\rm over}\geq P\f$. This function
   * should be simple enough that it does not depend on the evolution scale.
   * @param z   The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   */
  virtual double overestimateP(const double z, const IdList & ids) const = 0; 

  /**
   * Purely virtual method which should return
   * the ratio of the splitting function to the overestimate, i.e.
   * \f$P(z,\tilde{q}^2)/P_{\rm over}(z)\f$.
   * @param z   The energy fraction.
   * @param t   The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param mass Whether or not to include the mass dependent terms
   */
  virtual double ratioP(const double z, const Energy2 t, const IdList & ids,
			const bool mass, const RhoDMatrix & rho) const = 0;

  /**
   * Purely virtual method which should return the indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$.
   * @param z         The energy fraction.
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   *                  
   */
  virtual double integOverP(const double z, const IdList & ids, 
			    unsigned int PDFfactor=0) const = 0; 

  /**
   * Purely virtual method which should return the inverse of the 
   * indefinite integral of the 
   * overestimated splitting function, \f$P_{\rm over}\f$ which is used to
   * generate the value of \f$z\f$.
   * @param r Value of the splitting function to be inverted
   * @param ids The PDG codes for the particles in the splitting.
   * @param PDFfactor Which additional factor to include for the PDF
   *                  0 is no additional factor,
   *                  1 is \f$1/z\f$, 2 is \f$1/(1-z)\f$ and 3 is \f$1/z/(1-z)\f$
   */
  virtual double invIntegOverP(const double r, const IdList & ids, 
			       unsigned int PDFfactor=0) const = 0;
  //@}

  /**
   * Method to calculate the azimuthal angle for forward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param The azimuthal angle, \f$\phi\f$.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiForward(const double z, const Energy2 t, const IdList & ids,
		     const RhoDMatrix &) = 0;

  /**
   * Method to calculate the azimuthal angle for backward evolution
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @return The weight
   */
  virtual vector<pair<int,Complex> > 
  generatePhiBackward(const double z, const Energy2 t, const IdList & ids,
		      const RhoDMatrix &) = 0;

  /**
   * Calculate the matrix element for the splitting
   * @param z The energy fraction
   * @param t The scale \f$t=2p_j\cdot p_k\f$.
   * @param ids The PDG codes for the particles in the splitting.
   * @param phi The azimuthal angle, \f$\phi\f$.
   * @param timeLike Whether timelike or spacelike, affects inclusive of mass terms
   */
  virtual DecayMEPtr matrixElement(const double z, const Energy2 t, 
				   const IdList & ids, const double phi,
                                   bool timeLike) = 0;

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
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  Energy calculateScale(double z, Energy pt, IdList ids,unsigned int iopt);

  /**
   *  Scale choice
   */
  bool pTScale() const {
    return scaleChoice_ == 2 ? angularOrdered() : scaleChoice_ == 0;
  }

  /**
   * Method which should make the proper colour connection 
   * between the emitting parent and the branching products.
   * @param parent The parent for the branching
   * @param first  The first  branching product
   * @param second The second branching product
   * @param partnerType The type of evolution partner
   * @param back Whether this is foward or backward evolution.
   */
  void colourConnection(tShowerParticlePtr parent, tShowerParticlePtr first,
			tShowerParticlePtr second, ShowerPartnerType partnerType, 
			const bool back) const;

  /**
   *  Functions to state scales after branching happens
   */
  //@{
  /**
   *  Sort out scales for final-state emission
   */
  void evaluateFinalStateScales(ShowerPartnerType type,
				Energy scale, double z,
				tShowerParticlePtr parent,
				tShowerParticlePtr first,
				tShowerParticlePtr second);
  /**
   *  Sort out scales for initial-state emission
   */
  void evaluateInitialStateScales(ShowerPartnerType type,
				  Energy scale, double z,
				  tShowerParticlePtr parent,
				  tShowerParticlePtr first,
				  tShowerParticlePtr second);

  /**
   *  Sort out scales for decay emission
   */
  void evaluateDecayScales(ShowerPartnerType type,
			   Energy scale, double z,
			   tShowerParticlePtr parent,
			   tShowerParticlePtr first,
			   tShowerParticlePtr second);
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
   *  Methods to provide the next value of the scale before the vetos
   *  are applied.
   */
  //@{
  /**
   *  Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @param enhance The radiation enhancement factor
   * @return False if scale less than minimum, true otherwise
   */
  bool guessTimeLike(Energy2 &t, Energy2 tmin, double enhance, double detune);

  /**
   * Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmax The maximum scale
   * @param minmass The minimum mass of the particle after the branching
   * @param enhance The radiation enhancement factor
   */
  bool guessDecay(Energy2 &t, Energy2 tmax,Energy minmass,
		  double enhance, double detune);

  /**
   * Value of the energy fraction and scale for space-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @param x Fraction of the beam momentum.
   * @param enhance The radiation enhancement factor
   */
  bool guessSpaceLike(Energy2 &t, Energy2 tmin, const double x,
		      double enhance, double detune);
  //@}

  /**
   *  Initialize the values of the cut-offs and scales
   * @param tmin The minimum scale
   * @param ids  The ids of the partics in the branching
   */
  void initialize(const IdList & ids,Energy2 &tmin);

  /**
   *  Phase Space veto member to implement the \f$\Theta\f$ function as a veto
   *  so that the emission is within the allowed phase space.
   * @param t  The scale
   * @param maxQ2 The maximum virtuality
   * @return true if vetoed
   */
  bool PSVeto(const Energy2 t);

  /**
   * Compute the limits on \f$z\f$ for time-like branching
   * @param scale The scale of the particle
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeTimeLikeLimits(Energy2 & scale);

  /**
   * Compute the limits on \f$z\f$ for space-like branching
   * @param scale The scale of the particle
   * @param x The energy fraction of the parton
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeSpaceLikeLimits(Energy2 & scale, double x);

protected:

  /**
   *  Methods to implement the veto algorithm to generate the scale of 
   *  the next branching
   */
  //@{
  /**
   * Value of the energy fraction and value of the scale for the veto algorithm
   * @param iopt The option for calculating z
   * @param ids The PDG codes of the particles in the splitting
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   * @param t1 The starting valoe of the scale
   * @param enhance The radiation enhancement factor
   * @param identical Whether or not the outgoing particles are identical
   * @param t_main rerurns the value of the energy fraction for the veto algorithm
   * @param z_main returns the value of the scale for the veto algorithm
   */
  virtual void guesstz(Energy2 t1,unsigned int iopt, const IdList &ids,
		       double enhance,bool ident,
		       double detune, Energy2 &t_main, double &z_main);

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
    return ratioP(z_, t, ids,mass,rho)/detune;
  }
  //@}

protected:

  /**
   *  Return the colour factor
   */
  double colourFactor() const;

  /**
   *  The limits of \f$z\f$ in the splitting
   */
  pair<double,double> zLimits() const {return zlimits_;};

public:

  /**
   *  Calculate the virtual masses for a branchings
   */
  const vector<Energy> & virtualMasses(const IdList & ids) {
  	return cutoff_->virtualMasses(ids);
  }

  /**
   *  The minimum pT2
   */
  Energy2 pT2min() { return cutoff_->pT2min(); }

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Sudakov1to2FormFactor & operator=(const Sudakov1to2FormFactor &) = delete;

private:

  /**
   *  Pointer to the coupling for this Sudakov form factor
   */
  SudakovCutOffPtr cutoff_;

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

private:
  
  /**
   *  The evolution scale, \f$\tilde{q}\f$.
   */
  Energy q_;

  /**
   *  The Ids of the particles in the current branching
   */
  IdList ids_;

  /**
   *  The masses of the particles in the current branching
   */
  vector<Energy> masses_;

  /**
   *  The mass squared of the particles in the current branching
   */
  vector<Energy2> masssquared_;

private:
  /**
   *  The choice of scale
   */
  unsigned int scaleChoice_;

  /**
   *   Enforce strict AO
   */
  bool strictAO_;

  /**
   *  The colour factor
   */
  double colourFactor_;

};

}

#endif /* HERWIG_Sudakov1to2FormFactor_H */
