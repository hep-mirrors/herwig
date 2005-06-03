// -*- C++ -*-
#ifndef HERWIG_EtaPiPiGammaDecayer_H
#define HERWIG_EtaPiPiGammaDecayer_H
// This is the declaration of the EtaPiPiGammaDecayer class.

#include "Herwig++/Utilities/Kinematics.h"
#include "Herwig++/Decay/DecayIntegrator.h"
#include "Herwig++/Decay/DecayPhaseSpaceMode.h"
#include "Herwig++/Utilities/GaussianIntegral.h"
#include "Herwig++/Utilities/Interpolator.h"
#include "EtaPiPiGammaDecayer.fh"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The <code>EtaPiPiGammaDecayer</code> class implements the decay of
 * the \f$\eta\f$ or \f$\eta'\f$ to \f$\pi^+\pi^-\gamma\f$ using either 
 * a VMD type model or a model using either the theoretical or experimental 
 * form of the Omnes function taken from hep-ph/0112150.
 *
 * The matrix element is given by
 * \f[\mathcal{M} = B(s_{+-},s_{+\gamma},s_{-\gamma})\epsilon^{\mu\nu\alpha\beta}\epsilon^*_{\mu}p_{+\nu}p_{-\alpha}p_{\gamma\beta}\f]
 * where \f$p_{+,-}\f$ are the momenta of the positively and negatively charged pions,
 * \f$p_{\gamma}\f$ is the momentum of the photon and \f$s_{ij} = (p_i+p_j)^2\f$.
 *
 *  The different models take
 *
 *  \f[B(s_{+-},s_{+\gamma},s_{-\gamma}) = 
 *  B_0\left(1+\frac32\frac{s_{+-}}{M^2_\rho-s_{+-}-iM_\rho\Gamma_\rho(s_{+-})}\right)\f]
 *  where \f$M_\rho\f$ and \f$\Gamma_\rho\f$ are the mass and running width 
 *  of the \f$\rho\f$
 *  respectively for the VMD model.
 *
 *  For the Omnes function case we take
 *
 *  \f[B(s_{+-},s_{+\gamma},s_{-\gamma}) = 
 *  B_0\left(1-c+c\frac{1+as_{+-}}{D_1(s_{+-})}\right)\f]
 *  either the experimental or analytic form of the Omnes function \f$D_1(s_{+-})\f$
 *  taken from hep-ph/0112150 can be used.
 *
 *  The coefficient \f$B_0\f$ is given in hep-ph/0112150. We use the values from this
 *  paper and use their default choice \f$c=1\f$, \f$a=\frac1{2M_\rho}\f$.
 *
 * @see DecayIntegrator
 * 
 */
class EtaPiPiGammaDecayer: public DecayIntegrator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  EtaPiPiGammaDecayer();

  /**
   * Copy-constructor.
   */
  inline EtaPiPiGammaDecayer(const EtaPiPiGammaDecayer &);

  /**
   * Destructor.
   */
  virtual ~EtaPiPiGammaDecayer();
  //@}

public:

  /**
   * Accept member which is called at initialization to see if this Decayer can
   * handle a given decay mode. This version checks the particles against the 
   * list of allowed incoming  and outgoing mesons.
   * @param dm The DecayMode
   * @return Whether the mode can be handled.
   */
  virtual bool accept(const DecayMode & dm) const;
  
  /**
   * For a given decay mode and a given particle instance, perform the
   * decay and return the decay products. This version uses PDG codes to
   * work which mode is being simulated and the generate member of the 
   * DecayIntegrator class for the phase-space  generation.
   * @param dm The DecayMode
   * @param part The Particle instant being decayed.
   * @return The vector of particles produced in the decay.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & part) const;

  /**
   * Return the matrix element squared for a given mode and phase-space channel.
   * @param vertex Output the information on the vertex for spin correlations
   * @param ichan The channel we are calculating the matrix element for. 
   * @param part The decaying Particle.
   * @param decay The particles produced in the decay.
   * @return The matrix element squared for the phase-space configuration.
   */
  double me2(bool vertex, const int ichan,const Particle & part,
	     const ParticleVector & decay) const;

  /**
   * Method to return an object to calculate the 3 body partial width.
   * @param dm The DecayMode
   * @return A pointer to a WidthCalculatorBase object capable of calculating the width
   */
  virtual WidthCalculatorBasePtr threeBodyMEIntegrator(const DecayMode & dm) const;

  /**
   * The matrix element to be integrated for the three-body decays as a function
   * of the invariant masses of pairs of the outgoing particles.
   * @param imode The mode for which the matrix element is needed.
   * @param q2 The scale, \e i.e. the mass squared of the decaying particle.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element
   */
  virtual double threeBodyMatrixElement(int imode,Energy2 q2, Energy2 s3,Energy2 s2,
					Energy2 s1,Energy m1,Energy m2,Energy m3);

  /**
   * Output the setup information for the particle database
   */
  void dataBaseOutput(ofstream &);

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
  static ClassDescription<EtaPiPiGammaDecayer> initEtaPiPiGammaDecayer;

  /**
   * Private and non-existent assignment operator.
   */
  EtaPiPiGammaDecayer & operator=(const EtaPiPiGammaDecayer &);

private:

  /**
   * The analytic Omnes function, \f$D_1^{\rm anal}(s)\f$.
   * @param s The scale \f$s\f$.
   * @return The analytic Omnes function.
   */
  inline Complex analyticOmnes(Energy2 s) const;

  /**
   * The experimental Omnes function, \f$D_1^{\rm exp}(s)\f$.
   * @param s The scale \f$s\f$.
   * @return The experimental Omnes function.
   */
  inline Complex experimentalOmnes(Energy2 s) const;

private:

  /**
   * the pion decay constant, \f$F_\pi\f$.
   */
  Energy _fpi;

  /**
   * the PDG code for the incoming particle
   */
  vector<int> _incoming;

  /**
   * Coupling for the decay, \f$B_0\f$.
   */
  vector<double> _coupling;

  /**
   * The maximum weight
   */
  vector<double> _maxweight;

  /**
   * The option for the energy dependence of the prefactor
   */
  vector<int> _option;

  /**
   * The constants for the omnes function form.
   */
  InvEnergy2 _aconst;

  /**
   * The constants for the Omnes function form.
   */
  double _cconst;

  /**
   * The \f$\rho\f$ mass and width .
   */
  Energy _mrho,_rhowidth;

  /**
   * Constant for the running \f$rho\f$ width.
   */
  double _rhoconst;

  /**
   * The \f$m_\pi\f$.
   */
  Energy _mpi;

  /**
   * Use local values of the parameters.
   */
  bool _localparameters;

  /**
   *  Energy values for the experimental data on the phase shift
   */
  vector<Energy> _energy;

  /**
   *  Experimental values of the phase shift
   */
  vector<double> _phase;

  /**
   *  Energy values for the interpolation table for the Omnes function.
   */
  vector<Energy> _Omnesenergy;

  /**
   * Real part of the Omnes function for the interpolation table
   */
  vector<InvEnergy2> _Omnesfunctionreal;

  /**
   * Imaginary part of the Omnes function for the interpolation table
   */
  vector<InvEnergy2> _Omnesfunctionimag;

  /**
   * set up of the interpolation table
   */
  bool _initialize;

  /**
   * Number of points for the intepolation of the experimental Omnes function
   */
  unsigned int _npoints;

  /**
   * Interpolators for the experimental Omnes function.
   */
  mutable Interpolator *_Oreal,*_Oimag;

  /**
   *  Cut-off parameter for the integral of the experimental function
   */ 
  Energy _epscut;

  /**
   *  size parameters for the output
   */
  unsigned int _nsize[2];
 };

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

template <>
/**
 * The following template specialization informs ThePEG about the
 * base class of EtaPiPiGammaDecayer.
 */
struct BaseClassTrait<Herwig::EtaPiPiGammaDecayer,1> {
    /** Typedef of the base class of EtaPiPiGammaDecayer. */
  typedef Herwig::DecayIntegrator NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::EtaPiPiGammaDecayer>
  : public ClassTraitsBase<Herwig::EtaPiPiGammaDecayer> {
  /** Return the class name.*/
  static string className() { return "Herwig++::EtaPiPiGammaDecayer"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwSMDecay.so"; }

};

}



namespace Herwig {
using namespace Genfun;

  /** \ingroup Decay
   *
   *  A simple class used by the EtaPiPiGammaDecayer for the integration of the 
   * experimental Omnes function.
   *
   * @see EtaPiPiGammaDecayer
   *
   */
class OmnesFunction: public Genfun:: AbsFunction {

public:		   

  /**
   *  Function composition
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const; 

  /**
   * clone
   */
  OmnesFunction *clone() const; 

private:                               

  /**
   * clone
   */
  virtual AbsFunction *_clone() const;
    
public:
  
  /**
   * constructor with data as vectors
   */
  OmnesFunction(Interpolator *,Energy2);

  /**
   * Copy constructor
   */
  OmnesFunction(const OmnesFunction &right);
  
  /**
   * Standard ctors and dtor.
   */
  virtual ~OmnesFunction();
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double argument) const;

  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}

  /**
   *  set the scale
   */
  void setScale(Energy2);

private:
  
  /**
   * It is illegal to assign a function
   */
  const OmnesFunction & operator=(const OmnesFunction &right);
  
private:
  
  /**
   *  The interpolator
   */
  Interpolator *_interpolator;

  /**
   *  The scale and precision.
   */
  Energy2 _s,_precision; 
};
  

}

#include "EtaPiPiGammaDecayer.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "EtaPiPiGammaDecayer.tcc"
#endif

#endif /* HERWIG_EtaPiPiGammaDecayer_H */
