// -*- C++ -*-
#ifndef HERWIG_BaryonThreeQuarkModelFormFactor_H
#define HERWIG_BaryonThreeQuarkModelFormFactor_H
//
// This is the declaration of the BaryonThreeQuarkModelFormFactor class.
//

#include "BaryonFormFactor.h"
#include "BaryonThreeQuarkModelFormFactor.fh"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The BaryonThreeQuarkModelFormFactor class implements the 
   *  form factors for the semi-leptonic decay of baryons containing a heavy quark
   *  from the relativistic three-quark model calculation of PRD56, 348.
   *
   *  As the only formulae in the paper are for the heavy-to-heavy i.e. bottom
   *  to charm decay this there are the only modes included, although the paper
   *  also includes charm decays and bottom decays to light quarks.
   *
   *  The form factors are calculated by numerical computing the integrals from
   *  PRD56, 348 to obtaqin the coefficients for the expansion of the form factors.
   *
   * @see BaryonFormFactor
   */

class BaryonThreeQuarkModelFormFactor: public BaryonFormFactor {

  /**
   *  friend class for the integration of the expansion coefficents
   */
  friend class BaryonCFunction;

  /**
   *  friend class for the integration to give the partial width
   */
  friend class BaryonCMatrixElement;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline BaryonThreeQuarkModelFormFactor();

  /**
   * Copy constructor
   */
  inline BaryonThreeQuarkModelFormFactor(const BaryonThreeQuarkModelFormFactor &);

  /**
   * Destructor
   */
  virtual ~BaryonThreeQuarkModelFormFactor();
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

public:

  /** @name Form Factors */
  //@{
  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac12\f$ baryon.  
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming baryon.
   * @param id1 The PDG code of the outgoing baryon.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param f1v The form factor \f$F^V_1\f$.
   * @param f2v The form factor \f$F^V_2\f$.
   * @param f3v The form factor \f$F^V_3\f$.
   * @param f1a The form factor \f$F^A_1\f$.
   * @param f2a The form factor \f$F^A_2\f$.
   * @param f3a The form factor \f$F^A_3\f$.
   */
  virtual void SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
					  Energy m0, Energy m1,
					  Complex & f1v,Complex & f2v,Complex & f3v,
					  Complex & f1a,Complex & f2a,Complex & f3a);

  /**
   * The form factor for the weak decay of a spin \f$\frac12\f$ baryon to a 
   * spin \f$\frac32\f$ baryon.  
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form factor list.
   * @param id0 The PDG code of the incoming baryon.
   * @param id1 The PDG code of the outgoing baryon.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param g1v The form factor \f$G^V_1\f$.
   * @param g2v The form factor \f$G^V_2\f$.
   * @param g3v The form factor \f$G^V_3\f$.
   * @param g4v The form factor \f$G^V_4\f$.
   * @param g1a The form factor \f$G^A_1\f$.
   * @param g2a The form factor \f$G^A_2\f$.
   * @param g3a The form factor \f$G^A_3\f$.
   * @param g4a The form factor \f$G^A_4\f$.
   */
  virtual void SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
					       Energy m0, Energy m1,
					       Complex & g1v,Complex & g2v,Complex & g3v,
					       Complex & g4v,Complex & g1a,Complex & g2a,
					       Complex & g3a,Complex & g4a);
  //@}

  /**
   *  Output the information required to reproduce the object for the particle
   *  properties database
   */
  virtual void dataBaseOutput(ofstream&);

protected:

  /** @name Function needed to calculate the form factors */
  //@{
  /**
   * Returns the function \f$\Phi_N\f$ function of PRD56, 348 as a function
   * of \f$\omega\f$. 
   */
  inline vector<double> phiFunction(double);

  /**
   * The integral of a power of the the \f$S\f$ function of PRD56, 348 with respect
   * to \f$\theta\f$. This is used to 
   * calculate the coefficients of the expansion to compute the form factors.
   * The integral is calculated by computing a low power of the integrand and then
   * using recursion relations to calculate the pwoer requested.
   * @param y The \f$y\f$ variable of the function which is integrate over numerically.
   * @param N The power to which the integrand is raised.
   * @param SNm2 The integral with the function raised to the power \f$N/2-1\f$.
   * @param SN The integral with the function raised to the power \f$N\f$.  
   */
  inline void SN(double y,int N,double & SNm2,double & SN);
  // the integral of the S function

  /**
   *  The integrand for the coefficients of the expansion. This is a function of the
   * integration variable \f$x\f$ which is chosen to transform the integrand over \f$y\f$
   * which is from \f$0\f$ to \f$\infty\f$ to an integral between 0 and 1. This means
   * that \f$y=\frac{1-x}{x}\f$.
   * @param x The integration variable.
   */
  inline double integrandC(double x);
  // the integrand for the coefficients

  //@}

  /**
   *  The integrand for the semi-analytic calculation of the semi-leptonic width.
   *  This is mainly included for testing purposes.
   * @param omega The \f$\omega\f$ parameter of the heavy quark form-factors.
   * @param m0 The mass of the incoming baryon.
   * @param m1 The mass of the outgoing baryon.
   * @param type The type of the decay 
   * @param imass The baryonic mass parameter to use.
   */
  inline double widthIntegrand(double omega,Energy m0, Energy m1, int type, int imass);

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
   * Describe an abstract base class with persistent data.
   */
  static ClassDescription<BaryonThreeQuarkModelFormFactor> 
  initBaryonThreeQuarkModelFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  BaryonThreeQuarkModelFormFactor & operator=(const BaryonThreeQuarkModelFormFactor &);

private:

  /** @name Parameters for the form factors */
  //@{
  /**
   * Initialization of the expansion coefficients
   */
  bool _initialize;

  /**
   * Order of the expansion for the form factors.
   */
  unsigned int _order;

  /**
   *  Mass of the light quarks used in the calculation of the form factors.
   */
  Energy _mlight;

  /**
   *  Mass of the strange quark used in the calculation of the form factors.
   */
  Energy _mstrange;

  /**
   *  The heavy quark \f$\Lambda_Q\f$ parameter for the calculation of the form factors.
   */
  Energy _LambdaQ;

  /**
   *  The \f$\Lambda_{qq}\f$ parameter for the calculation of the form factors.
   */
  Energy _Lambdaqq;

  /**
   *  The \f$\Lambda_{sq}\f$ parameter for the calculation of the form factors.
   */
  Energy _Lambdasq;


  /**
   *  The \f$\Lambda_{ss}\f$ parameter for the calculation of the form factors.
   */
  Energy _Lambdass;

  /**
   *  Coefficients for the expansion of the \f$\xi(\omega)\f$ form factor.
   */
  vector<double> _C0;

  /**
   *  Coefficients for the expansion of the \f$\xi_1(\omega)\f$ form factor.
   */
  vector<double> _C1;

  /**
   *  Coefficients for the expansion of the \f$\xi_2(\omega)\f$ form factor.
   */
  vector<double> _C2;

  /**
   *  Coefficient of the first term in the integrand for the coefficient calculation.
   */
  double _a;

  /**
   *  Coefficient of the second term in the integrand for the coefficient calculation.
   */
  double _b;

  /**
   *  The \f$\mu^2\f$ parameter for the coefficient calculation.
   */
  double _mu2;

  /**
   *  the order of the coefficient being calculated.
   */
  int _N;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * BaryonThreeQuarkModelFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::BaryonThreeQuarkModelFormFactor,1> {
  /** Typedef of the base class of BaryonThreeQuarkModelFormFactor. */
   typedef Herwig::BaryonFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BaryonThreeQuarkModelFormFactor class.
 */
template <>
 struct ClassTraits<Herwig::BaryonThreeQuarkModelFormFactor>
  : public ClassTraitsBase<Herwig::BaryonThreeQuarkModelFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig++::BaryonThreeQuarkModelFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwFormFactor.so"; }
};

}

// class for the integration of the coefficients
namespace Herwig {
using namespace Genfun;
using namespace ThePEG;

  /** \ingroup Decay
   *  This is a function using the CLHEP Genfun class whiches can access the integrandC
   * member of the BaryonThreeQuarkModelFormFactor class. This function can then
   * be integrated to give the coefficients.
   */
class BaryonCFunction : public Genfun::AbsFunction {

public:		   

  /**
   *  Function composition
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const; 

  /**
   * clone
   */
  BaryonCFunction *clone() const; 

private:                               

  /**
   * clone
   */
  virtual AbsFunction *_clone() const;

public:
 
/** @name Standard constructors and destructors. */
/**
 *  The constructor
 */
 BaryonCFunction(BaryonThreeQuarkModelFormFactorPtr);
  
  /**
   *  The destructor
   */
  virtual ~BaryonCFunction();
  
  /**
   * The copy constructor
   */
  BaryonCFunction(const BaryonCFunction &right);
  //@}

  /**
   *  Retreive the function value
   */
  virtual double operator ()(double argument) const;

  /**
   *  Retreive the function value
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
private:
  
  /**
   * Non-existant assignment operator. It is illegal to assign a function
   */
  const BaryonCFunction & operator=(const BaryonCFunction &right);
  
private:
  
  /**
   *  A pointer to the form factor to supply the integrand.
   */
  BaryonThreeQuarkModelFormFactorPtr _formFactor;

};
}


// class for the integration of the coefficients
namespace Herwig {
using namespace Genfun;
using namespace ThePEG;

/** \ingroup Decay
 *  This function function uses the CLHEP Genfun classes which can access the 
 *  widthIntegrand member of the BaryonThreeQuarkModelFormFactor class. This function
 * can then be integrated to give the partial width.
 *
 */
class BaryonCMatrixElement : public Genfun::AbsFunction {

public:		   

  /**
   *  Function composition
   */
  virtual FunctionComposition operator()(const AbsFunction &function) const; 

  /**
   * clone
   */
  BaryonCMatrixElement *clone() const; 

private:                               

  /**
   * clone
   */
  virtual AbsFunction *_clone() const;

public:

/**
 * The constructor
 */
 BaryonCMatrixElement(BaryonThreeQuarkModelFormFactorPtr,Energy,Energy,int,int);

/** @name Standard constructors and destructors. */
  /**
   * The destructor
   */
  virtual ~BaryonCMatrixElement();

  /**
   *  The copy constructor 
   */
  BaryonCMatrixElement(const BaryonCMatrixElement &right);
  //@}

  /**
   *  Retreive the function value
   */
  virtual double operator ()(double argument) const;

  /**
   *  Retreive the function value
   */
  virtual double operator ()(const Argument & a) const {return operator() (a[0]);}
  
private:
  
  /**
   * Non-existant assignment operator. It is illegal to assign a function
   */
  const BaryonCMatrixElement & operator=(const BaryonCMatrixElement &right);
  
private:
  
  /**
   *  A pointer to the form factor to supply the integrand.
   */
  BaryonThreeQuarkModelFormFactorPtr _formFactor;

  /** The mass of the incoming baryon. */
  Energy _m0;

  /** The mass of the outgoing baryon. */
  Energy _m1;

  /** The type of decay being integrated. */
  int _type;

  /** The baryonic mass parameter. */
  int _mass;

};
}


#include "BaryonThreeQuarkModelFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "BaryonThreeQuarkModelFormFactor.tcc"
#endif

#endif /* HERWIG_BaryonThreeQuarkModelFormFactor_H */
