// -*- C++ -*-
#ifndef HERWIG_TwoMesonRhoKStarCurrent_H
#define HERWIG_TwoMesonRhoKStarCurrent_H
// This is the declaration of the TwoMesonRhoKStarCurrent class.

#include "WeakDecayCurrent.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "TwoMesonRhoKStarCurrent.fh"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/**  \ingroup Decay
 *
 *  Weak current for the production of two mesons via the \f$\rho\f$ or \f$K^*\f$
 *  resonances.
 *  These currents are taken from tau decays.
 *
 *  The current takes the form
 *
 *  \f[J^\mu = \frac{\sqrt{2}}{\sum_k\alpha_k}\left((p_1-p_2)^\mu-\frac{(p_1-p_2)\cdot q}{q^2}q^\mu))\right)
 *   \sum_k \alpha_k B_{R_k}(q^2)
 *  \f]
 *  where
 *  - \f$p_{1,2}\f$ are the momenta of the outgoing mesons,
 *  - \f$q=p_1+p_2\f$,
 *  - \f$B_{R_k}(q^2)\f$ is the Breit-Wigner distribution for the intermediate vector
 *    meson \f$R_k\f$.
 *  - \f$\alpha_k\f$ is the weight for the resonance.
 *
 *   The Breit-Wigner term is summed over the \f$\rho\f$ or \f$K^*\f$ resonances that
 *   can contribute to a given decay.
 *
 *  The models of either Kuhn and Santamaria (Z. Phys. C48, 445 (1990))
 *  or Gounaris and Sakurai Phys. Rev. Lett. 21, 244 (1968) are supported for the
 *  shape of the Breit-Wigner distribution. The mixing parameters
 *  are taken from Phys.Rev.D61:112002,2000 (CLEO) for the decay \f$\pi^\pm\pi^0\f$ and
 *  the CLEO version of TAUOLA for the \f$K\pi\f$ decays.
 *
 * @see WeakDecayCurrent.
 * 
 *  \author Peter Richardson
 *
 */
class TwoMesonRhoKStarCurrent: public WeakDecayCurrent {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline TwoMesonRhoKStarCurrent();

  /**
   * Copy constructor
   */
  inline TwoMesonRhoKStarCurrent(const TwoMesonRhoKStarCurrent &);

  /**
   * Destructor
   */
  virtual ~TwoMesonRhoKStarCurrent();
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

public:

  /** @name Methods for the construction of the phase space integrator. */
  //@{ 
  /**
   * Complete the construction of the decay mode for integration.
   * This version just adds the intermediate resonances, two outgoing mesons
   * and photon.
   * @param icharge The total charge of the outgoing particles in the current.
   * @param imode   The mode in the current being asked for.
   * @param mode    The phase space mode for the integration
   * @param iloc    The location of the of the first particle from the current in
   *                the list of outgoing particles.
   * @param ires    The location of the first intermediate for the current.
   * @param phase   The prototype phase space channel for the integration.
   * @param upp     The maximum possible mass the particles in the current are
   *                allowed to have.
   * @return Whether the current was sucessfully constructed.
   */
  virtual bool createMode(int icharge,unsigned int imode,DecayPhaseSpaceModePtr mode,
			  unsigned int iloc,unsigned int ires,
			  DecayPhaseSpaceChannelPtr phase,Energy upp);

  /**
   * The particles produced by the current. This just returns the two pseudoscalar
   * mesons and the photon.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version returns the hadronic current described above.
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector>  current(bool vertex, const int imode,
						     const int ichan,Energy & scale,  
						     const ParticleVector & decay) const;

  /**
   * Accept the decay. Checks the particles are the allowed mode.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * @param id The id's of the particles in the current.
   * @return The number of the mode
   */
  virtual unsigned int decayMode(vector<int> id);

  /**
   * Output the information for the database
   */
  virtual void dataBaseOutput(ofstream &) const;

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
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
   * @throws RebindException if no cloned object was found for a given pointer.
   */
  inline virtual void rebind(const TranslationMap & trans)
    throw(RebindException);

  /**
   * Return a vector of all pointers to Interfaced objects used in
   * this object.
   * @return a vector of pointers.
   */
  inline virtual IVector getReferences();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<TwoMesonRhoKStarCurrent> initTwoMesonRhoKStarCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  TwoMesonRhoKStarCurrent & operator=(const TwoMesonRhoKStarCurrent &);

private:

  /** @name Breit-Wigner functions */
  //@{
  /**
   * \f$p\f$-wave breit wigner for form-factors
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param imodel Which of the two models for the Breit-Wigner shape to use.
   * @param itype  Whether this is a \f$\rho\f$ or \f$K^*\f$ resonance.
   * @param ires   Which of the different multiplets to use.
   * @return The value of the Breit-Wigner distribution.
   */
  inline Complex BreitWigner(Energy2 q2, unsigned int imodel,unsigned int itype,
			     unsigned int ires) const;
  
  /**
   * The \f$d\f$ parameter in the GS model of the propagator.
   * @param ires Which of the \f$\rho\f$ or \f$K^*\f$ resonances to use
   * @return The \f$d\f$ parameter
   */
  inline double GSModelDParameter(const unsigned int ires)const ;
  
  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the GS model for the propagator
   * evaluated at \f$q^2=m^2\f$.
   * @param ires Which of the \f$\rho\f$ or \f$K^*\f$ resonances to use
   * @return The value of \f$\frac{dh}{dq^2}\f$ at \f$q^2=m^2\f$.
   */
  inline double GSModeldhdq2Parameter(const unsigned int ires)const ;
  
  /**
   * The \f$h\f$ function in the GS model.
   * @param ires Which of the \f$\rho\f$ or \f$K^*\f$ resonances to use
   * @param q The scale \f$q\f$.
   * @return The value of the \f$h\f$ function at scale \f$q\f$.
   */
  inline double GSModelhFunction(const unsigned int ires, const Energy q)const ;
  
  /**
   * The momentum of the decay products for the running width.
   * @param ires Which of the \f$\rho\f$ or \f$K^*\f$ resonances to use
   * @param q The scale \f$q\f$.
   * @return The momenta of the decay products in the rest frame of the
   *         intermediate resonance.
   */
  inline Energy pcm(const unsigned int ires, const Energy q) const;
  //@}

private:

  /**
   * Weights for the different \f$\rho\f$ resonances in the current, \f$\alpha_k\f$.
   */
  vector<double> _piwgt;

  /**
   * Weights for the different \f$K^*\f$ resonances in the current, \f$\alpha_k\f$.
   */
  vector<double> _kwgt;

  /**
   * Model to use for the \f$\rho\f$ propagator.
   */
  int _pimodel;

  /**
   * Model to use for the \f$K^*\f$ propagator.
   */
  int _Kmodel;

  /**
   * Option not to use the physical masses and widths for the \f$\rho\f$.
   */
  bool _rhoparameters;

  /**
   * The masses of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhomasses;

  /**
   * The widths of the \f$\rho\f$ resonances.
   */
  vector<Energy> _rhowidths;
  
  /**
   * Option not to use the physical masses and widths for the \f$K^*\f$.
   */
  bool _Kstarparameters;

  /**
   *  The masses of the \f$K^*\f$ resonances.
   */
  vector<Energy> _Kstarmasses;

  /**
   *  The masses of the \f$K^*\f$ resonances.
   */
  vector<Energy> _Kstarwidths;

  /**
   * Parameters for the Breit-Wigners
   */
  //@{
  /**
   * The masses of the resonances
   */
  vector<Energy> _mass;

  /**
   * The widths of the resonances
   */
  vector<Energy> _width;

  /**
   * squares of the masses
   */
  vector<Energy2> _mass2;

  /**
   * product of the mass and the width
   */
  vector<Energy2> _massw;

  /**
   * Masses of the decay products for the momentum calculation.
   */
  vector<Energy> _massa,_massb;

  /**
   * The decay for the on-shell mass momentum
   */
  vector<Energy>  _mom;

  /**
   * The function \f$\frac{dh}{dq^2}\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<InvEnergy2> _dhdq2;

  /**
   * The function \f$h\f$ at \f$q^2=m^2\f$ for the GS form of the
   *  Breit-Wigner
   */
  vector<InvEnergy2> _hm2;

  /**
   * The \f$d\f$ parameter  for the GS form of the
   *  Breit-Wigner
   */
  vector<double> _dparam;
  //@}
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of TwoMesonRhoKStarCurrent.
 */
template <>
struct BaseClassTrait<Herwig::TwoMesonRhoKStarCurrent,1> {
  /** Typedef of the base class of VectorMesonCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::TwoMesonRhoKStarCurrent>
  : public ClassTraitsBase<Herwig::TwoMesonRhoKStarCurrent> {
  /** Return the class name. */
  static string className() { return "Herwig++::TwoMesonRhoKStarCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwWeakCurrent.so"; }

};

}

#include "TwoMesonRhoKStarCurrent.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TwoMesonRhoKStarCurrent.tcc"
#endif

#endif /* HERWIG_TwoMesonRhoKStarCurrent_H */
