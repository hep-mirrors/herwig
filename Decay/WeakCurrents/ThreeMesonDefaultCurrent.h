// -*- C++ -*-
#ifndef HERWIG_ThreeMesonDefaultCurrent_H
#define HERWIG_ThreeMesonDefaultCurrent_H
//
// This is the declaration of the ThreeMesonDefaultCurrent class.
//
#include "ThreeMesonCurrentBase.h"
#include "ThreeMesonDefaultCurrent.fh"
#include "Herwig++/Utilities/NewInterpolator.h"
#include "Herwig++/Decay/ThreeBodyIntegrator.h"
#include "Herwig++/Utilities/Kinematics.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 * The ThreeMesonDefaultCurrent class implements the currents from Z.Phys.C58:445 (1992),
 * this paper uses the form from Z.Phys.C48:445 (1990) for the \f$a_1\f$ width and
 * is the default model in TAUOLA.
 *
 *  The following three meson modes are implemented.
 *
 * - \f$    \pi^-  \pi^-    \pi^+ \f$, (imode=0)
 * - \f$    \pi^0  \pi^0    \pi^- \f$, (imode=1)
 * - \f$    K^-   \pi^-    K^+ \f$, (imode=2)
 * - \f$    K^0   \pi^-    \bar{K}^0\f$, (imode=3)
 * - \f$    K^-   \pi^0    K^0 \f$, (imode=4)
 * - \f$    \pi^0  \pi^0    K^- \f$, (imode=5)
 * - \f$    K^-   \pi^-    \pi^+ \f$, (imode=6)
 * - \f$    \pi^-  \bar{K}^0  \pi^0 \f$, (imode=7)
 * - \f$    \pi^-  \pi^0    \eta \f$, (imode=8)
 *
 *  using the currents from TAUOLA
 *
 *
 * @see ThreeMesonCurrentBase,
 * @see WeakDecayCurrent.
 * @see Defaulta1MatrixElement
 * 
 */
class ThreeMesonDefaultCurrent: public ThreeMesonCurrentBase {

  /**
   * The matrix element for the running \f$a_1\f$ width is a friend to 
   * keep some members private.
   */
  friend class Defaulta1MatrixElement;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  ThreeMesonDefaultCurrent();

  /**
   * Copy constructor
   */
  inline ThreeMesonDefaultCurrent(const ThreeMesonDefaultCurrent &);
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
   * This version addes the mesons for the current
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
   * The particles produced by the current. This returns the mesons for the mode.
   * @param icharge The total charge of the particles in the current.
   * @param imode The mode for which the particles are being requested
   * @param iq The PDG code for the quark
   * @param ia The PDG code for the antiquark
   * @return The external particles for the current.
   */
  virtual PDVector particles(int icharge, unsigned int imode, int iq, int ia);
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;
  
protected:

  /**
   * the matrix element for the \f$a_1\f$ decay to calculate the running width
   * @param q2 The mass of the decaying off-shell \f$a_1\f$, \f$q^2\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param m1 The mass of the first  outgoing particle.
   * @param m2 The mass of the second outgoing particle.
   * @param m3 The mass of the third  outgoing particle.
   * @return The matrix element squared summed over spins.
   */
  inline double a1MatrixElement(Energy2 q2, Energy2 s3,Energy2 s2,Energy2 s1,
				Energy m1,Energy m2,Energy m3);

  /**
   * Can the current handle a particular set of mesons. 
   * As this current includes all the allowed modes this is always true.
   */
  virtual bool acceptMode(int) const;

  /**
   * Calculate the form factor for the current. Implements the form factors
   * described above.
   * @param ichan The phase space channel
   * @param imode The mode
   * @param q2 The scale \f$q^2\f$ for the current.
   * @param s1 The invariant mass squared of particles 2 and 3, \f$s_1=m^2_{23}\f$.
   * @param s2 The invariant mass squared of particles 1 and 3, \f$s_2=m^2_{13}\f$.
   * @param s3 The invariant mass squared of particles 1 and 2, \f$s_3=m^2_{12}\f$.
   * @param F1 The form factor \f$F_1\f$.
   * @param F2 The form factor \f$F_2\f$.
   * @param F3 The form factor \f$F_3\f$.
   * @param F4 The form factor \f$F_4\f$.
   * @param F5 The form factor \f$F_5\f$.
   */
  virtual void calculateFormFactors(const int ichan,const int imode,
				    Energy2 q2,Energy2 s1,Energy2 s2,Energy2 s3,
				    Complex&F1,Complex&F2,Complex&F3,
				    Complex&F4,Complex&F5) const;

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
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object to the begining of the run phase.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ThreeMesonDefaultCurrent> initThreeMesonDefaultCurrent;

  /**
   * Private and non-existent assignment operator.
   */
  ThreeMesonDefaultCurrent & operator=(const ThreeMesonDefaultCurrent &);

private:
  
  /**
   * The \f$\rho\f$ Breit-Wigner for the \f$F_{1,2,3}\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  inline Complex BrhoF123(Energy2 q2,int ires) const;

  /**
   * The \f$\rho\f$ Breit-Wigner for the \f$F_5\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  inline Complex BrhoF5(Energy2 q2,int ires) const;

  /**
   * The \f$K^*\f$ Breit-Wigner for the \f$F_{1,2,3}\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  inline Complex BKstarF123(Energy2 q2,int ires) const;

  /**
   * The \f$K^*\f$ Breit-Wigner for the \f$F_5\f$ form factors.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @param ires Which \f$\rho\f$ multiplet
   * @return The Breit-Wigner 
   */
  inline Complex BKstarF5(Energy2 q2,int ires) const;
  
  /**
   * Mixed Breit Wigner for the \f$F_5\f$ form factor
   * @param s1 The scale \f$s_1\f$.
   * @param s2 The scale \f$s_2\f$.
   * @param ires Which resonances to use
   * @return The mixed Breit-Wigner
   */
  inline Complex FKrho(Energy2 s1,Energy2 s2,int ires) const;
  
  /**
   * \f$a_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  inline Complex a1BreitWigner(Energy2 q2) const;
  
  /**
   * The \f$K_1\f$ Breit-Wigner
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The Breit-Wigner
   */
  inline Complex K1BreitWigner(Energy2 q2) const;
  
  /**
   * The \f$a_1\f$ running width
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner
   * @return The \f$a_1\f$ running width.
   */
  inline Energy a1width(Energy2 q2) const ;
  
  /**
   * Initialize the \f$a_1\f$ running width
   * @param iopt Initialization option (-1 full calculation, 0 set up the interpolation)
   */
  inline void inita1width(int iopt);

  /**
   * Breit-Wigners for the \f$\rho\f$ and \f$K^*\f$.
   * @param q2 The scale \f$q^2\f$ for the Breit-Wigner.
   * @param itype The type of Breit-Wigner, \e i.e. which masses and widths to use.x
   * @param ires Which multiplet to use.
   */
  inline Complex rhoKBreitWigner(Energy2 q2,unsigned int itype,unsigned int ires) const;

private:
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_{1,2,3}\f$ form factors.
   */
  vector<double> _rhoF123wgts;

  /**
   * Parameters for the \f$K^*\f$ Breit-Wigner in the
   * \f$F_{1,2,3}\f$ form factors.
   */
  vector<double> _KstarF123wgts;
  
  /**
   * Parameters for the \f$\rho\f$ Breit-Wigner in the
   * \f$F_5\f$ form factors.
   */
  vector<double> _rhoF5wgts;

  /**
   * Parameters for the \f$K^*\f$ Breit-Wigner in the
   * \f$F_5\f$ form factors.
   */
  vector<double> _KstarF5wgts;
  
  /**
   * The relative weight of the \f$\rho\f$ and \f$K^*\f$ where needed.
   */
  double _rhoKstarwgt;
  
  /**
   * The \f$a_1\f$ width for the running \f$a_1\f$ width calculation.
   */
  vector<Energy>  _a1runwidth;

  /**
   * The \f$q^2\f$ for the running \f$a_1\f$  width calculation.
   */
  vector<Energy2> _a1runq2;


  /**
   * The interpolator for the running \f$a_1\f$ width calculation.
   */
  NewInterpolatorPtr _a1runinter;

  /**
   * Initialize the running \f$a_1\f$ width.
   */
  bool _initializea1;
  
  /**
   * The mass of the \f$a_1\f$ resonances.
   */
  Energy _a1mass;

  /**
   * The width of the \f$a_1\f$ resonances.
   */
  Energy _a1width;

  /**
   * The mass of the \f$aK1\f$ resonances.
   */
  Energy _K1mass;

  /**
   * The width of the \f$K_1\f$ resonances.
   */
  Energy _K1width;

  /**
   * The pion decay constant, \f$f_\pi\f$.
   */
  Energy _fpi;

  /**
   * The pion mass
   */
  Energy _mpi;

  /**
   * The kaon mass
   */
  Energy _mK;

  /**
   * use local values of the \f$\rho\f$ masses and widths
   */
  bool _rhoparameters;

  /**
   * The \f$\rho\f$ masses for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _rhoF123masses;

  /**
   * The \f$\rho\f$ masses for the \f$F_5\f$ form factors.
   */
  vector<Energy> _rhoF5masses;

  /**
   * The \f$\rho\f$ widths for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _rhoF123widths;

  /**
   * The \f$\rho\f$ widths for the \f$F_5\f$ form factors.
   */
  vector<Energy> _rhoF5widths;
  
  /**
   * use local values of the \f$K^*\f$ resonances masses and widths
   */
  bool _Kstarparameters;

  /**
   * The \f$K^*\f$ masses for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _KstarF123masses;

  /**
   * The \f$K^*\f$ masses for the \f$F_5\f$ form factors.
   */
  vector<Energy> _KstarF5masses;

  /**
   * The \f$K^*\f$ widths for the \f$F_{1,2,3}\f$ form factors.
   */
  vector<Energy> _KstarF123widths;

  /**
   * The \f$K^*\f$ widths for the \f$F_5\f$ form factors.
   */
  vector<Energy> _KstarF5widths;
  
  /**
   * Use local values of the \f$a_1\f$ parameters
   */
  bool _a1parameters;
  
  /**
   * Use local values of the \f$K_1\f$ parameters
   */
  bool _K1parameters;
  
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * The following template specialization informs ThePEG about the
 * base class of ThreeMesonDefaultCurrent.
 */
template <>
 struct BaseClassTrait<Herwig::ThreeMesonDefaultCurrent,1> {
  /** Typedef of the base class of ThreeMesonDefaultCurrent. */
  typedef Herwig::ThreeMesonCurrentBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<Herwig::ThreeMesonDefaultCurrent>
  : public ClassTraitsBase<Herwig::ThreeMesonDefaultCurrent> {
  /** Return the class name. */
  static string className() { return "Herwig++::ThreeMesonDefaultCurrent"; }
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwWeakCurrents.so"; }

};

}

#include "CLHEP/GenericFunctions/AbsFunction.hh"
namespace Herwig {
using namespace Genfun;
using namespace ThePEG; 

/** \ingroup Decay
 *
 * Definitions of the functions to be integrated to give the running
 * function to return the matrix element for the \f$a_1\f$ decay to be
 * integrated to give the \f$a_1\f$ running width
 *
 * @see ThreeMesonDefaultCurrent
 *
 */
class Defaulta1MatrixElement : public Genfun::AbsFunction {
        
public:
  
  /**
   * FunctionComposition operator
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const;
  
  /**
   * Clone method
   */
   Defaulta1MatrixElement *clone() const;

private:

  /**
   * Clone method
   */
  virtual AbsFunction *_clone() const;
    
public:

  /**
   * Constructor
   */
  Defaulta1MatrixElement(ThreeMesonDefaultCurrentPtr);

  /**
   *  The number of variables, in thsi case 7
   */  
  virtual unsigned int dimensionality() const ;     
  
  /**
   * Copy constructor
   */
  Defaulta1MatrixElement(const Defaulta1MatrixElement &right);
  
  /**
   * Retreive function value
   */
  virtual double operator ()(double) const {return 0.;}
  
  /**
   * Retreive function value
   */
  virtual double operator ()(const Argument & a) const ;
  
  /**
   *   set the scale 
   */
  inline void setQ2(Energy2);
  
  
private:
  
  /**
   * It is illegal to assign a function
   */
  const Defaulta1MatrixElement & 
  operator=(const Defaulta1MatrixElement &right);
  
private:

  /**
   * The current
   */
  ThreeMesonDefaultCurrentPtr _decayer;
};
}

#include "ThreeMesonDefaultCurrent.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ThreeMesonDefaultCurrent.tcc"
#endif

#endif /* THEPEG_ThreeMesonDefaultCurrent_H */
