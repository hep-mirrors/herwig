// -*- C++ -*-
#ifndef HERWIG_SudakovFormFactor_H
#define HERWIG_SudakovFormFactor_H
//
// This is the declaration of the SudakovFormFactor class.
//

#include "Herwig++/Shower/ShowerConfig.h"
#include "ThePEG/Interface/Interfaced.h"
#include "Herwig++/Utilities/GlobalParameters.h"
#include "SplittingFunction.h"
#include "Herwig++/Shower/ShowerAlpha.h"
#include "Herwig++/Shower/ShowerIndex.h"
#include "Herwig++/Shower/ShowerVariables.h"
#include "SudakovFormFactor.fh"

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
 *  given the previous value $\tilde{q}_i$
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
 *  - tVeto()           member to ensure the scale is above the minimum value.
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
 */
class SudakovFormFactor: public Interfaced {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline SudakovFormFactor();

  /**
   * The copy constructor.
   */
  inline SudakovFormFactor(const SudakovFormFactor &);

  /**
   * The destructor.
   */
  virtual ~SudakovFormFactor();
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
   * @param revOrd is used (when it is not equal to the default, false, value)
   *        only for final-state branching of a decaying on-shell particle.
   */
  virtual Energy generateNextTimeBranching(const Energy startingScale,
				           const IdList &ids,
				           const bool revOrd = false);

  /**
   * Return the scale of the next space-like branching. If there is no 
   * branching then it returns Energy().
   * @param startingScale starting scale for the evolution
   * @param ids The PDG codes of the particles in the splitting
   * @param beam    Pointer to the particleData for the 
   *                beam
   * @param pdf Pointer to the PDF object
   * @param x The fraction of the beam momentum
   * @param revOrd is used (when it is not equal to the default, false, value)
   *        only for final-state branching of a decaying on-shell particle.
   */
  virtual Energy generateNextSpaceBranching(const Energy startingScale,
		                            const IdList &ids,
					    const tcPDPtr beam,
					    const tcPDFPtr pdf,
					    double x,
					    const bool revOrd = false);


  //@}

  /**
   * This virtual method is defined as empty, and it should be
   * overriden only for those derived Sudakov form factor classes
   * that use lookup tables for numerical evaluations, rather
   * than using the Monte Carlo rejection (veto) method.
   * This method is called once, during initialization, by
   * the SplittingGenerator. 
   * General methods, usable for any type of Sudakov form factor
   * that override this method, should be provided in this
   * class in the protected session.
   */
  virtual void setupLookupTables();

  /**
   *  Methods to provide public access to the private member functions
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
  inline ShowerIndex::InteractionType interactionType();

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
   *  Methods to provide the next value of the scale before the vetos
   *  are applied.
   */
  //@{
  /**
   *  Value of the scale for the veto algorithm
   * @param t1 The starting valoe of the scale
   * @param z0 Lower limit on \f$z\f$.
   * @param z1 Upper limit on \f$z\f$.
   * @param IS Whether this is for the space-like (true) or time-like (false)
   * shower.
   */
  Energy2 guesst (Energy2 t1,double z0,double z1,bool IS=false) const;

  /**
   * Value of the energy fraction for the veto algorithm
   * @param z0 Lower limit on \f$z\f$.
   * @param z1 Upper limit on \f$z\f$.
   */
  double guessz (double z0, double z1) const;

  /**
   *  Value of the energy fraction and scale for time-like branching
   * @param z The energy fraction
   * @param z0 Lower limit on \f$z\f$.
   * @param z1 Upper limit on \f$z\f$.
   * @param t  The scale
   * @param t0 The cut-off on the scale
   * @param kinCutoff The kinematic cut-off
   * @param glueEmits Whether or not the emitting particle is a gluon.
   */
  void guessTimeLike(double &z, double &z0, double &z1,
		     Energy2 &t, const Energy2 t0,
		     const Energy kinCutoff, const bool glueEmits) const;

  /**
   * Value of the energy fraction and scale for space-like branching
   * @param z The energy fraction
   * @param z0 Lower limit on \f$z\f$.
   * @param z1 Upper limit on \f$z\f$.
   * @param t  The scale
   * @param t0 The cut-off on the scale
   * @param kinCutoff The kinematic cut-off
   * @param x Fraction of the beam momentum.
   */
  void guessSpaceLike(double &z, double &z0, double &z1,
	     Energy2 &t, const Energy2 t0,
	     const Energy kinCutoff, const double x) const;
  //@}

  /**
   *  Initialize the values of the cut-offs
   * @param t0 The cut-off scale
   * @param tmin The minimum scale
   * @param tmax The maximum scale
   * @param kinCutoff The kinematic cut-off
   * @param m The mass of the emitting particle
   */
  void initialize(Energy2 &t0, Energy2 &tmin, Energy2 tmax, 
		  Energy &kinCutoff, Energy m);

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
  bool alphaSVeto(const Energy2 pt2) const;

  /**
   *  Veto on the PDF for the initial-state shower
   * @param z The energy fraction
   * @param t The scale
   * @param x The fraction of the beam momentum
   * @param pdf Pointer to the PDF object
   * @param parton0 Pointer to the particleData for the 
   *                new parent (this is the particle we evolved back to)
   * @param parton1 Pointer to the particleData for the 
   *                original particle
   * @param beam    Pointer to the particleData for the 
   *                beam
   */
  bool PDFVeto(const double z, const Energy2 t, const double x,
	       const tcPDFPtr pdf, const tcPDPtr parton0, 
	       const tcPDPtr parton1, const tcPDPtr beam) const;

  /**
   *  Phase Space veto
   * @param z The energy fraction
   * @param z0 Lower limit on \f$z\f$.
   * @param z1 Upper limit on \f$z\f$.
   * @param t  The scale
   * @param t0 The cut-off value of the scale
   * @param tmin The minimum value of the scale
   * @param kinCutoff The kinematic cut-off
   * @param glueEmits Whether or not the emitting particle is a gluon.
   * @return true if vetoed
   */
  bool PSVeto(const double z, const double z0, const double z1,
	      const Energy2 t, const Energy2 tmin, const Energy2 t0,
	      const Energy kinCutoff, bool glueEmits);

  /**
   *  The veto on the splitting function.
   * @param z The energy fraction
   * @param t The scale
   * @param ids The PDG codes of the particles in the splitting 
   * @return true if vetoed
   */
  bool SplittingFnVeto(const double z, const Energy2 t, 
		       const IdList &ids) const;

  /**
   * Veto on the scale
   * @param t Scale which is reset to -1 if fails veto
   * @param tmin The minimum value of the scale 
   * @return true if vetoed
   */
  bool tVeto(Energy2 &t, const Energy2 tmin) const;
  //@}

  /**
   * Toy model: returns \f$q\f$ according to powerlike (\f$q^p\f$) distribution
   * with cutoff \f$q_0\f$, \f$q_{min} < q_0 < q_{max}\f$, 
   * where  \f$q_{min}\f$ is chosen such that the 
   * probability for a first branching is 1-R. 
   * \f$z\f$ is selected according to distributions 
   * below depending on the choice of znorm.
   * @param znorm Distribution to use 
   *    \f$1/(1-z)\f$ with \f$z_0=m/q\f$ and \f$z_0 < z < 1-z_0\f$ (true) 
   * or \f$z^2+(1-z)^2\f$ with \f$z_0 < z < 1\f$ (false)
   * @param p The power
   * @param R Parameter controlling the probability of the first branching,
   *          which is \f$1-R\f$.
   * @param q0 The cut-off scale
   * @param qmax The maximum value of the scale
   * @param q The scale.
   * @param z The energy fraction.
   */
  void get_qz (bool znorm, double p, double R, Energy q0, Energy qmax, 
	       Energy &q, double &z);
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
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

protected:

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
  double _z;

  /**
   *  The azimuthal angle
   */
  double _phi;
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
  static string className() {return "Herwig++::SudakovFormFactor"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the SudakovFormFactor class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "SudakovFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SudakovFormFactor.tcc"
#endif

#endif /* HERWIG_SudakovFormFactor_H */
