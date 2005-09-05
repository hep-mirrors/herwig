// -*- C++ -*-
#ifndef HERWIG_FourPionNovosibirskCurrent_H
#define HERWIG_FourPionNovosibirskCurrent_H
//
// This is the declaration of the FourPionNovosibirskCurrent class.
//
#include "WeakDecayCurrent.h"
#include "FourPionNovosibirskCurrent.fh"
#include "Herwig++/Utilities/Interpolator.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>FourPionNovosibirskCurrent</code> class implements the decay of the weak 
 * current to 4 pions using the hadronic currents of
 * Comput. Phys. Commun. 146: 139-153, 2002,
 * which is a model based on the \f$e^+e^-\to4\pi\f$ data from Novosibirsk.
 *
 * It should be noted that there were a large number of mistakes in this paper which 
 * were corrected in hep-ph/0312240.
 *
 * @see WeakDecayCurrent
 * @see FourPionDefaultMatrixElement
 * 
 * \author Peter Richardson
 * 
 */
class FourPionNovosibirskCurrent: public WeakDecayCurrent {

  /**
   * The FourPionDefaultMatrixElement class is a friend so it can perform the
   * integration.
   */
  friend class FourPionDefaultMatrixElement;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  FourPionNovosibirskCurrent();

  /**
   * Copy constructor
   */
  inline FourPionNovosibirskCurrent(const FourPionNovosibirskCurrent &);

  /**
   * Destructor
   */
  virtual ~FourPionNovosibirskCurrent();
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
   * This version constructs the four pion current.
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
   * The particles produced by the current. This returns the four pions for the
   * current.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Hadronic current. This version calculates the four pion current described above.
   * @param vertex Construct the information needed for spin correlations
   * @param imode The mode
   * @param ichan The phase-space channel the current is needed for.
   * @param scale The invariant mass of the particles in the current.
   * @param decay The decay products
   * @return The current. 
   */
  virtual vector<LorentzPolarizationVector>  current(bool vertex, const int imode,
						     const int ichan,Energy & scale, 
						     const ParticleVector & decay) const;

  /**
   * Accept the decay. Checks this is one of the four pion modes.
   * @param id The id's of the particles in the current.
   * @return Can this current have the external particles specified.
   */
  virtual bool accept(vector<int> id);

  /**
   * Return the decay mode number for a given set of particles in the current. 
   * Works out which four pion mode this is.
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
  inline virtual void doinit() throw(InitException);

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
  static ClassDescription<FourPionNovosibirskCurrent> initFourPionNovosibirskCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  FourPionNovosibirskCurrent & operator=(const FourPionNovosibirskCurrent &);
      
protected:
  
  /**
   * The matrix element to evaluate the \f$a_1\f$ running width.
   * @param q2 The mass of the decaying off-shell \f$a_1\f$, \f$q^2\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element squared summed over spins.
   */
  inline double a1MatrixElement(Energy2 q2, Energy2 s3,Energy2 s2, Energy2 s1,
				Energy m1,Energy m2,Energy m3);

  /**
   * Initialize the \f$a_1\f$ width.
   * @param iopt Initialization option
   *  (-1 is full initialization and 0 sets up the interpolator for the running width)
   */
  inline void inita1width(int iopt);
  
  /**
   * Form foactor for the \f$a_1\f$ vertex.
   * @param q2 The scale \f$q^2\f$.
   * @return The \f$a_1\f$ form factor.
   */
  inline double a1FormFactor(Energy2 q2) const;

  /**
   * Breit-Wigner for the \f$\sigma\f$ meson
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$\sigma\f$ meson
   */
  inline Complex sigmaBreitWigner(Energy2 q2) const;

  /**
   * The \f$a_1\f$ breit wigner.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$a_1\f$.
   */
  inline Complex a1BreitWigner(Energy2 q2) const;

  /**
   * The Breit-Wigner for the \f$\omega\f$.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$\omega\f$.
   */
  inline Complex omegaBreitWigner(Energy2 q2) const;

  /**
   * The Breit-Wigner for the \f$\rho\f$.
   * @param q2 The scale \f$q^2\f$.
   * @return The Breit-Wigner for the \f$\rho\f$.
   */
  inline Complex rhoBreitWigner(Energy2 q2) const;

  /**
   * Return the \f$a_1\f$ running width.
   * @param q2 The scale \f$q^2\f$.
   * @return The running width.
   */
  inline Energy a1width(Energy2 q2) const ;

  /**
   * The \f$t_1\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @return The current \f$t_1\f$.
   */
  inline LorentzPolarizationVector t1(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
				      Lorentz5Momentum & q3,Lorentz5Momentum & q4) const;

  /**
   * The \f$t_2\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @return The current \f$t_2\f$.
   */
  inline LorentzPolarizationVector t2(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
				      Lorentz5Momentum & q3,Lorentz5Momentum & q4) const;

  /**
   * The \f$t_3\f$ current used in calculating the current.
   * @param q1 The first momentum.
   * @param q2 The first momentum.
   * @param q3 The first momentum.
   * @param q4 The first momentum.
   * @return The current \f$t_3\f$.
   */
  inline LorentzPolarizationVector t3(Lorentz5Momentum & q1,Lorentz5Momentum & q2,
				      Lorentz5Momentum & q3,Lorentz5Momentum & q4) const;

  /**
   * The G functions of hep-ph/0201149
   * @param q2 The scale \f$q^2\f$.
   * @param ichan Which of the four pion channels this is for.
   * @return The G function.
   */
  inline double gFunction(Energy2 q2, int ichan) const;

  /**
   * The d parameter in \f$\rho\f$ the propagator.
   */
  inline Energy2 DParameter() const ;

  /**
   * The \f$\frac{dh}{dq^2}\f$ function in the rho propagator evaluated at \f$q^2=m^2\f$.
   */
  inline double dhdq2Parameter() const ;
  
  /**
   * The h function in the \f$\rho\f$ propagator.
   * @param q The scale.
   * @return The h function.
   */
  inline Energy2 hFunction(const Energy q)const ;

private:
  
  /**
   * Interpolating functions for the G functions of hep-ph/0201149
   */
  //@{
  /**
   * The interpolator for the \f$\omega\f$ current.
   */
  Interpolator *_Fomega;

  /**
   * The interpolator for the three charged pion \f$a_1\f$ current. 
   */
  Interpolator *_Fthreec;

  /**
   * The interpolator for the one   charged pion \f$a_1\f$ current.
   */
  Interpolator *_Fonec;

  /**
   * The interpolator for the \f$\sigma\f$ current.
   */
  Interpolator *_Fsigma;
  //@}

  /**
   * the pion mass
   */
  Energy _mpi;

  /**
   * The mass of the \f$\rho\f$ for the current.
   */
  Energy _rhomass;

  /**
   * The mass of the \f$a_1\f$ for the current.
   */
  Energy _a1mass;

  /**
   * The mass of the \f$\omega\f$ for the current.
   */
  Energy _omegamass;

  /**
   * The mass of the \f$\sigma\f$ for the current.
   */
  Energy _sigmamass;

  /**
   * The width for the \f$\rho\f$.
   */
  Energy _rhowidth;

  /**
   *  The \f$a_1\f$ width
   */
  Energy _a1width;

  /**
   *  The \f$\omega\f$ width.
   */
  Energy _omegawidth;

  /**
   *  The \f$\sigma\f$ width.
   */
  Energy _sigmawidth;

  /**
   * Mass for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intmass;

  /**
   * Width for the intermediate in the phase-space, this is a technical parameter to
   * improve the phase-space integration efficiency.
   */
  Energy _intwidth; 

  /**
   * The \f$z\f$ \f$\sigma\f$ coupling.
   */
  Complex _zsigma;

  /**
   * The magnitude of the \f$z\f$ \f$\sigma\f$ coupling.
   */
  double _zmag;

  /**
   * The phase of the \f$z\f$ \f$\sigma\f$ coupling.
   */
  double _zphase;

  /**
   * The mass parameter for the \f$a_1\f$ form-factor.
   */
  Energy2 _lambda2;

  /**
   * The inverse of the mass parameter for the \f$a_1\f$ form-factor.
   */
  InvEnergy2 _onedlam2;

  /**
   *  The physical \f$a_1\f$ mass divided by the mass parameter in the 
   *  \f$a_1\f$ form-factor.
   */
  double _a1massolam2;

  /**
   * The momentum of the  pions in on-shell \f$\sigma\f$ decay which is used
   * in the calculation of the running \f$\sigma\f$ width.
   */
  Energy _psigma;

  /**
   * The pion mass squared.
   */
  Energy2 _mpi2;

  /**
   *  The h function evaluated at the \f$\rho\f$ mass.
   */
  Energy2 _hm2;

  /**
   * The d parameter for the \f$\rho\f$ width.
   */
  Energy2 _rhoD;

  /**
   * The momentum of the pions produced in on-shell \f$rho\f$ decay.
   */
  Energy _prho;

  /**
   * \f$\frac{dh}{dq^2}\f$ evaluates at \f$q^2=m^2\f$ for the \f$\rho\f$.
   */
  double _dhdq2m2;

  /**
   * Magic number for the omega current.
   */
  InvEnergy2 _aomega;

  /**
   * Magic number for the three charged pion current.
   */
  InvEnergy2 _athreec;

  /**
   * Magic number for the one charged pion current
   */
  InvEnergy2 _aonec;

  /**
   * Magic number for the omega current.
   */
  double _bomega;

  /**
   * Magic number for the three charged pion current.
   */
  double _bthreec;

  /**
   * Magic number for the one charged pion current
   */
  double _bonec;

  /**
   * Magic number for the omega current.
   */
  Energy _comega;

  /**
   * Magic number for the three charged pion current.
   */
  Energy _cthreec;

  /**
   * Magic number for the one charged pion current
   */
  Energy _conec;

  /**
   * magic numbers for the running omega width
   */
  vector<double> _omegaparam;

  /**
   * whether or not to initialize the calculation of the \f$a_1\f$ width
   */
  bool _initializea1;

  /**
   * use local values of the particle masses
   */
  bool _localparameters;

  /**
   * The widths for the interpolation table for the running \f$a_1\f$ width.
   */
  vector<Energy> _a1runwidth;

  /**
   * The \f$q^2\f$ values for the interpolation table for the running \f$a_1\f$ width.
   */
  vector<Energy2> _a1runq2;

  /**
   * The interpolator for the running \f$a_1\f$ width.
   */
  Interpolator *_a1runinter;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of FourPionNovosibirskCurrent.
 */
template <>
 struct BaseClassTrait<Herwig::FourPionNovosibirskCurrent,1> {
   /** Typedef of the base class of FourPionNovosibirskCurrent. */
  typedef Herwig::WeakDecayCurrent NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::FourPionNovosibirskCurrent>
  : public ClassTraitsBase<Herwig::FourPionNovosibirskCurrent> {
  /** Return the class name.*/
  static string className() { return "Herwig++::FourPionNovosibirskCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwWeakCurrent.so"; }
  
};

}

#include "CLHEP/GenericFunctions/AbsFunction.hh"
namespace Herwig {

using namespace Genfun;
using namespace ThePEG; 

/** \ingroup Decay
 *
 * Definitions of the functions to be integrated to give the running \f$a_1\f$ width
 * function to return the matrix element for the \f$a_1\f$ decay to be
 * integrated to give the \f$a_1\f$ running width in the FourPionNovosibirskCurrent.
 *
 * @see FourPionNovosibirskCurrent
 *
 */
class FourPionDefaultMatrixElement : public Genfun::AbsFunction {
      
public:
  
  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
  FourPionDefaultMatrixElement *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;

public:

  /**
   * Constructor
   */
  FourPionDefaultMatrixElement(FourPionNovosibirskCurrentPtr);
  
  /**
   * Destructor
   */
  virtual ~FourPionDefaultMatrixElement();
  
  /**
   *  The number of variables in this case 7.
   */
  virtual unsigned int dimensionality() const ;     
  
  /**
   * Copy constructor
   */
  FourPionDefaultMatrixElement(const FourPionDefaultMatrixElement &right);
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double) const {return 0.;}

  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const ;

  /**
   * set the mass of the decay \f$a_1\f$.
   */  
  inline void setQ2(Energy2);
  
  
private:
  
  /**
   * It is illegal to assign a function
   */
  const FourPionDefaultMatrixElement & 
  operator=(const FourPionDefaultMatrixElement &right);
  
  /**
   * the decayer
   */
private:

  /**
   *  pointer to the current
   */
  Ptr<Herwig::FourPionNovosibirskCurrent>::pointer _decayer;
};

}


#include "FourPionNovosibirskCurrent.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FourPionNovosibirskCurrent.tcc"
#endif

#endif /* HERWIG_FourPionNovosibirskCurrent_H */

