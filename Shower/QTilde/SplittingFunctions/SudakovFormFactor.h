// -*- C++ -*-
//
// SudakovFormFactor.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
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
#include "Herwig/Shower/QTilde/ShowerConfig.h"
#include "Herwig/Shower/QTilde/Kinematics/ShowerKinematics.fh"
#include "ThePEG/EventRecord/RhoDMatrix.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig/Shower/ShowerAlpha.h"
#include "SudakovFormFactor.fh"

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 * Enum to define the possible types of colour structure which can occur in
 * the branching.
 */
enum ColourStructure {Undefined=0,
                      TripletTripletOctet  = 1, OctetOctetOctet       = 2,
                      OctetTripletTriplet  = 3, TripletOctetTriplet   = 4,
                      SextetSextetOctet    = 5, TripletTripletSinglet = 6,
		      OctetOctetSinglet    = 7, Epsilon               = 8,
		      OctetSinglet         = 9,
		      ChargedChargedNeutral=-1,
		      ChargedNeutralCharged=-2,
		      NeutralChargedCharged=-3,
		      EW=-4};

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

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SudakovFormFactor() : interactionType_(ShowerInteraction::UNDEFINED),
                        angularOrdered_(true),
                        colourStructure_(Undefined),
                        pdfMax_(35.0), pdfFactor_(0)
  {}
  //@}

public:

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
					      double enhance, double detuning) = 0;

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
						 double detuning) = 0;

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
						 double detuning) = 0;
  //@}

  /**
   *  Azimuthal angle generation
   */
  //@{
  /**
   * Generate the azimuthal angle of the branching for forward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiForward(ShowerParticle & particle,const IdList & ids,
				    ShoKinPtr kinematics, const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for backward evolution
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiBackward(ShowerParticle & particle,const IdList & ids,
				     ShoKinPtr kinematics, const RhoDMatrix & rho)=0;

  /**
   *  Generate the azimuthal angle of the branching for ISR in decays
   * @param particle The branching particle
   * @param ids The PDG codes of the particles in the branchings
   * @param The Shower kinematics
   */
  virtual double generatePhiDecay(ShowerParticle & particle,const IdList & ids,
				  ShoKinPtr kinematics, const RhoDMatrix & rho)=0;
  //@}
  
public:

  /**
   *  Return the colour structure
   */
  ColourStructure colourStructure() const {return colourStructure_;}
  
  /**
   *  Method to check the colours are correct
   */
  bool checkColours(const IdList & ids) const;

  /**
   * Return the pointer to the ShowerAlpha object.
   */
  tShowerAlphaPtr alpha() const { return alpha_; }

public:

  /**
   *  Methods to return the interaction type and order for the Sudakov form factor
   */
  //@{
  /**
   *  Purely virtual method which should determine whether this splitting
   *  function can be used for a given set of particles.
   *  @param ids The PDG codes for the particles in the splitting.
   */
  virtual bool accept(const IdList & ids) const = 0;

  /**
   *  Method to return the evolution scale given the
   *  transverse momentum, \f$p_T\f$ and \f$z\f$.
   */
  virtual Energy calculateScale(double z, Energy pt, IdList ids,unsigned int iopt) = 0;

  /**
   *  Return the type of the interaction
   */
  ShowerInteraction interactionType() const {return interactionType_;}

  /**
   *  Whether or not the interaction is angular ordered
   */
  bool angularOrdered() const {return angularOrdered_;}
  //@}

protected:

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

  /**
   *  The PDF factor
   */
  unsigned pdfFactor() const {return pdfFactor_;}

  /**
   * Maximum value of the PDF weight
   */
  double pdfMax() const { return pdfMax_;}
  
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
  bool PDFVeto(const Energy2 t, const double x, const double z,
	       const tcPDPtr parton0, const tcPDPtr parton1,
	       tcBeamPtr beam) const;
  /**
   * The PDF veto ratio
   */
  double PDFVetoRatio(const Energy2 t, const double x, const double z,
               const tcPDPtr parton0, const tcPDPtr parton1,
               tcBeamPtr beam,double factor) const;

  /**
   *   Set the PDF
   */
  void setPDF(tcPDFPtr pdf, Energy scale) {
    pdf_ = pdf;
    freeze_ = scale;
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
  virtual double alphaSVetoRatio(Energy2 pt2,double factor) const;

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

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovFormFactor & operator=(const SudakovFormFactor &) = delete;

private:

  /**
   *  The interaction type for the splitting function.
   */
  ShowerInteraction interactionType_;

  /**
   *  Whether or not this interaction is angular-ordered
   */
  bool angularOrdered_;

  /**
   *  The colour structure
   */
  ColourStructure colourStructure_;

  /**
   * List of the particles this Sudakov is used for to aid in setting up
   * interpolation tables if needed
   */
  vector<IdList> particles_;

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

  /**
   * Maximum value of the PDF weight
   */
  double pdfMax_;
  
  /**
   *  Option for the inclusion of a factor \f$1/(1-z)\f$ in the PDF estimate
   */
  unsigned pdfFactor_;
  //@}

  /**
   *  Pointer to the coupling for this Sudakov form factor
   */
  ShowerAlphaPtr alpha_;

};

}

#endif /* HERWIG_SudakovFormFactor_H */
