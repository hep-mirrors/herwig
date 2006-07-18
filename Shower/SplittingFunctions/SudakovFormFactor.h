// -*- C++ -*-
#ifndef HERWIG_SudakovFormFactor_H
#define HERWIG_SudakovFormFactor_H
//
// This is the declaration of the SudakovFormFactor class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "SplittingGenerator.fh"
#include "ThePEG/Repository/UseRandom.h"
#include "Herwig++/Shower/ShowerConfig.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "SplittingFunction.h"
#include "Herwig++/Shower/Couplings/ShowerIndex.h"
#include "Herwig++/Shower/ShowerVariables.h"
#include "SudakovFormFactor.fh"
#include <cassert>

namespace Herwig {

using namespace ThePEG;

/** \ingroup Shower
 *
 *  This is the definition of the Sudakov form factor class.
 *  The implementation of the form factors in this class uses the veto algorithm
 *  rather than look-up tables. This differs from FORTRAN HERWIG and the class is
 *  designed so that classes which inherit from this one can use look-up tables 
 *  if needed.
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
 *  - SplittingFnVeto() member for the veto on the value of the splitting function
 *  - PSVeto()          member to implement the \f$\Theta\f$ function as a veto
 *                      so that the emission is within the allowed phase space.
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
 *  @see SplittingFunction
 *  @see ShowerAlpha
 *  @see SplittingGenerator
 *  @see \ref SudakovFormFactorInterfaces "The interfaces"
 *  defined for SudakovFormFactor.
 */
class SudakovFormFactor: public Interfaced {

/**
 *  Splitting Generator is a friend so that it can set the ShowerVariables
 */
friend class SplittingGenerator;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SudakovFormFactor();
  //@}

public:

  /**
   *  Methods to provide public access to the private member variables
   */
  //@{
  /** 
   * Return the pointer to the SplitFun object.
   */
  inline tSplittingFnPtr splittingFn() const;

  /**
   * Return the pointer to the ShowerAlpha object.
   */
  inline tShowerAlphaPtr alpha() const;
  //@}

  /**
   * Methods to provide public access to the variables for the shower kinematics
   * which are kept internally in this class and are generated by a call to 
   * generateNextTimeBranching or generateNextSpaceBranching. These variables
   * cannot be set externally directly but only via a call to the 
   * generateNextTimeBranching or generateNextSpaceBranching methods. The first
   * action of these methods is to clear the values they got in the previous
   * call to the same method. The lifetime 
   * of the values of these kinematics variables is therefore between
   * to successive call to generateNextBranching.
   * Finally, at the moment these variables are meaninful
   * only for a 1->2 splitting, but in future other variables
   * could be added as well for describing also a 1->3 splitting.
   */
  //@{
  /**
   *  The type of interaction
   */
  inline ShowerIndex::InteractionType interactionType() const;

  /**
   *  The energy fraction
   */
  inline double z() const;

  /**
   *  The azimuthal angle
   */
  inline double phi() const;

  /**
   *  The evolution scale
   */
  inline Energy qtilde() const;

  /**
   *  The resolution scale
   */
  inline Energy resScale() const;

  /**
   *  The kinematic scale
   */
  inline Energy kinScale() const;

  /**
   *  The transverse momentum
   */
  inline Energy pT() const;
  //@}

public:

  /**
   *  Methods to get the scale of the next branching
   */
  //@{
  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   */
  virtual Energy generateNextTimeBranching(const Energy startingScale,
				           const IdList &ids);

  /**
   * Return the scale of the next time-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param stoppingScale stopping scale for the evolution
   * @param minmass The minimum mass allowed for the spake-like particle.
   * @param ids The PDG codes of the particles in the splitting
   */
  virtual Energy generateNextDecayBranching(const Energy startingScale,
					    const Energy stoppingScale,
					    const Energy minmass,
					    const IdList &ids);

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param x The fraction of the beam momentum
   */
  virtual Energy generateNextSpaceBranching(const Energy startingScale,
		                            const IdList &ids,double x);


  //@}

  /**
   * This virtual method is defined as empty, and it should be
   * overriden only for those derived Sudakov form factor classes
   * that use lookup tables for numerical evaluations, rather
   * than using the Monte Carlo rejection (veto) method.
   * This method is called once, during initialization, by
   * the SplittingGenerator. 
   * General methods, usable  any type of Sudakov form factor
   * that override this method, should be provided in this
   * class in the protected session.
   */
  virtual void setupLookupTables();

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
   *  Methods to provide the next value of the scale before the vetos
   *  are applied.
   */
  //@{
  /**
   *  Value of the scale for the veto algorithm
   * @param t1 The starting valoe of the scale
   * @param iopt The option for calculating t 
   * - 0 is final-state
   * - 1 is initial-state for the hard process
   * - 2 is initial-state for particle decays
   */
  inline Energy2 guesst (Energy2 t1,unsigned int iopt) const;

  /**
   * Value of the energy fraction for the veto algorithm
   */
  inline double guessz () const;

  /**
   *  Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @return False if scale less than minimum, true otherwise
   */
  bool guessTimeLike(Energy2 &t, Energy2 tmin) const;

  /**
   * Value of the energy fraction and scale for time-like branching
   * @param t  The scale
   * @param tmax The maximum scale
   * @param minmass The minimum mass of the particle after the branching
   */
  bool guessDecay(Energy2 &t, Energy2 tmax,Energy minmass) const;

  /**
   * Value of the energy fraction and scale for space-like branching
   * @param t  The scale
   * @param tmin The minimum scale
   * @param x Fraction of the beam momentum.
   */
  bool guessSpaceLike(Energy2 &t, Energy2 tmin, const double x) const;
  //@}

  /**
   *  Initialize the values of the cut-offs
   * @param tmin The minimum scale
   * @param tmax The maximum scale
   * @param m The mass of the emitting particle
   */
  inline void initialize(Energy2 &t0, Energy2 &tmin, Energy2 tmax, 
			 Energy &kinCutoff, Energy m);


  /**
   *  Initialize the values of the cut-offs and scales
   * @param tmin The minimum scale
   * @param ids  The ids of the partics in the branching
   */
  void initialize(const IdList & ids,Energy2 &tmin);

  /**
   * The various different vetos which need to be applied using the veto
   * algorithm 
   */
  //@{
  /**
   *  The veto on the coupling constant
   * @param pt2 The value of ther transverse momentum squared, \f$p_T^2\f$.
   * @return true if vetoed
   */
  inline bool alphaSVeto(const Energy2 pt2) const;

  /**
   *  Veto on the PDF for the initial-state shower
   * @param z The energy fraction
   * @param t The scale
   * @param x The fraction of the beam momentum
   * @param parton0 Pointer to the particleData for the 
   *                new parent (this is the particle we evolved back to)
   * @param parton1 Pointer to the particleData for the 
   *                original particle
   */
  bool PDFVeto(const double z, const Energy2 t, const double x,
	       const tcPDPtr parton0, const tcPDPtr parton1) const;

  /**
   *  Phase Space veto
   * @param t  The scale
   * @return true if vetoed
   */
  bool PSVeto(const Energy2 t);

  /**
   *  The veto on the splitting function.
   * @param z The energy fraction
   * @param t The scale
   * @param ids The PDG codes of the particles in the splitting 
   * @return true if vetoed
   */
  inline bool SplittingFnVeto(const double z, const Energy2 t, 
		       const IdList &ids) const;
  //@}

  /**
   * Compute the limits on \f$z\f$ for time-like branching
   * @param scale The scale of the particle
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeTimeLikeLimits(Energy2 & scale) const;

  /**
   * Compute the limits on \f$z\f$ for space-like branching
   * @param scale The scale of the particle
   * @param x The energy fraction of the parton
   * @return True if lower limit less than upper, otherwise false
   */
  bool computeSpaceLikeLimits(Energy2 & scale, double x) const;
  /**
   *  Set the ShowerVariables
   */
  inline void setShowerVariables(ShowerVarsPtr);
  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<SudakovFormFactor> initSudakovFormFactor;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SudakovFormFactor & operator=(const SudakovFormFactor &);

private:

  /**
   * Member variables to keep the shower kinematics information
   * generated by a call to generateNextTimeBranching or generateNextSpaceBranching
   */
  //@{
  /**
   *  The evolution scale, \f$\tilde{q}\f$.
   */
  Energy _q;

  /**
   *  The energy fraction
   */
  mutable double _z;

  /**
   *  The azimuthal angle
   */
  double _phi;

  /**
   *  The transverse momentum
   */
  Energy _pt;
  //@}

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
   *  Pointer to the shower variables
   */
  ShowerVarsPtr _variables;

  /**
   * Maximum value of the PDF weight
   */
  double _pdfmax;

  /**
   *  The Ids of the particles in the current branching
   */
  IdList _ids;

  /**
   *  The masses of the particles in the current branching
   */
  vector<Energy> _masses;

  /**
   *  The mass squared of the particles in the current branching
   */
  vector<Energy> _masssquared;

  /**
   *  Kinematic cut-off
   */
  Energy _kinCutoff;

  /**
   *  The limits of \f$z\f$ in the splitting
   */
 mutable pair<double,double> _zlimits;

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
