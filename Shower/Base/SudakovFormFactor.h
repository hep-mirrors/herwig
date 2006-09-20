// -*- C++ -*-
#ifndef HERWIG_SudakovFormFactor_H
#define HERWIG_SudakovFormFactor_H
//
// This is the declaration of the SudakovFormFactor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingFunction.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Couplings/ShowerIndex.h"
#include "Herwig++/Shower/SplittingFunctions/SplittingGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include <cassert>
#include "ShowerKinematics.h"
#include "SudakovFormFactor.fh"

namespace Herwig {

using namespace ThePEG;

/**  \ingroup Shower
 *
 *  This is the definition of the Sudakov form factor class. In general this
 *  is the base class for the implementation of Sudakov form factors in Herwig++.
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
  inline SudakovFormFactor();

  /**
   *  Members to generate the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * @param enhance The radiation enhancement factor
   * defined.
   */
  virtual ShoKinPtr generateNextTimeBranching(const Energy startingScale,
					      const IdList &ids,const bool cc,
					      double enhance)=0;

  /**
   * Return the scale of the next space-like decay branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param enhance The radiation enhancement factor
   */
  virtual ShoKinPtr generateNextDecayBranching(const Energy startingScale,
					       const Energy stoppingScale,
					       const Energy minmass,
					       const IdList &ids,
					       const bool cc,
					       double enhance)=0;

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   * @param cc Whether this is the charge conjugate of the branching
   * defined.
   * @param beam The beam particle
   * @param enhance THe radiation enhancement factor
   */
  virtual ShoKinPtr generateNextSpaceBranching(const Energy startingScale,
					       const IdList &ids,double x,
					       const bool cc,double enhance,
					       Ptr<BeamParticleData>::transient_const_pointer beam)=0;
  //@}

  /**
   *  Methods to provide public access to the private member variables
   */
  //@{
  /** 
   * Return the pointer to the SplittingFunction object.
   */
  inline tSplittingFnPtr splittingFn() const;

  /**
   * Return the pointer to the ShowerAlpha object.
   */
  inline tShowerAlphaPtr alpha() const;

  /**
   *  The type of interaction
   */
  inline ShowerIndex::InteractionType interactionType() const;
  //@}

public:

  /**
   *  Methods to access the kinematic variables for the branching
   */
  //@{
  /**
   *  The energy fraction
   * @param Whether this is thge charge conjugate of the branching defined in the
   * splitting function.
   */
  inline double z() const;

  /**
   *  The azimuthal angle
   */
  inline double phi() const;

  /**
   *  The transverse momentum
   */
  inline Energy pT() const;
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

  /**
   *  Methods to implement the veto algorithm to generate the scale of 
   *  the next branching
   */
  //@{
  /**
   * Value of the energy fraction for the veto algorithm
   */
  inline double guessz () const;

  /**
   *  Value of the scale for the veto algorithm
   * @param t1 The starting valoe of the scale
   * @param iopt The option for calculating t 
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   * @param enhance The radiation enhancement factor
   * @param identical Whether or not the outgoing particles are identical
   */
  inline Energy2 guesst (Energy2 t1,unsigned int iopt,
			 double enhance, bool identical) const;

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
	       Ptr<BeamParticleData>::transient_const_pointer beam) const;

  /**
   *  The veto on the splitting function.
   * @param t The scale
   * @param ids The PDG codes of the particles in the splitting
   * @param mass Whether or not to use the massive splitting functions 
   * @return true if vetoed
   */
  inline bool SplittingFnVeto(const Energy2 t, 
			      const IdList &ids, 
			      const bool mass) const;

  /**
   *  The veto on the coupling constant
   * @param pt2 The value of ther transverse momentum squared, \f$p_T^2\f$.
   * @return true if vetoed
   */
  inline bool alphaSVeto(const Energy2 pt2) const;
  //@}

  /**
   *  Methods to set the kinematic variables for the branching
   */
  //@{
  /**
   *  The energy fraction
   */
  inline void z(double);

  /**
   *  The azimuthal angle
   */
  inline void phi(double);

  /**
   *  The transverse momentum
   */
  inline void pT(Energy);
  //@}

  /**
   *  Set/Get the limits on the energy fraction for the splitting
   */
  //@{
  /**
   * Get the limits
   */
  inline pair<double,double> zLimits() const;

  /**
   * Set the limits
   */
  inline void zLimits(pair<double,double> );
  //@}

  /**
   *  Set the particles in the splittings
   */
  void addSplitting(const IdList &);

  /**
   *  Access the potential branchings
   */
  inline vector<IdList> particles() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<SudakovFormFactor> initSudakovFormFactor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovFormFactor & operator=(const SudakovFormFactor &);

private:

  /**
   *  Pointer to the splitting function for this Sudakov form factor
   */
  SplittingFnPtr _splittingFn;

  /**
   *  Pointer to the coupling for this Sudakov form factor
   */
  ShowerAlphaPtr _alpha;

  /**
   * Maximum value of the PDF weight
   */
  double _pdfmax;

  /**
   * List of the particles this Sudakov is used for to aid in setting up
   * interpolation tables if needed
   */
  vector<IdList> _particles;

private:

  /**
   * Member variables to keep the shower kinematics information
   * generated by a call to generateNextTimeBranching or generateNextSpaceBranching
   */
  //@{
  /**
   *  The energy fraction
   */
  double _z;

  /**
   *  The azimuthal angle
   */
  double _phi;

  /**
   *  The transverse momentum
   */
  Energy _pT;
  //@}

private:

  /**
   *  The limits of \f$z\f$ in the splitting
   */
  pair<double,double> _zlimits;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SudakovFormFactor. */
template <>
struct BaseClassTrait<Herwig::SudakovFormFactor,1> {
  /** Typedef of the first base class of SudakovFormFactor. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SudakovFormFactor class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::SudakovFormFactor>
  : public ClassTraitsBase<Herwig::SudakovFormFactor> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::SudakovFormFactor"; }
  /**
   * The name of a file containing the dynamic library where the class
   * SudakovFormFactor is implemented. It may also include several, space-separated,
   * libraries if the class SudakovFormFactor depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "SudakovFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SudakovFormFactor.tcc"
#endif

#endif /* HERWIG_SudakovFormFactor_H */
