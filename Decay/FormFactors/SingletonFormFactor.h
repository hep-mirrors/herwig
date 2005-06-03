// -*- C++ -*-
#ifndef HERWIG_SingletonFormFactor_H
#define HERWIG_SingletonFormFactor_H
//
// This is the declaration of the SingletonFormFactor class.
//

#include "BaryonFormFactor.h"
#include "ThePEG/PDT/ParticleData.h"
#include "SingletonFormFactor.fh"
#include "CLHEP/GenericFunctions/AbsFunction.hh"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig++/Utilities/GaussianIntegral.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The SingletonFormFactor class implements the form-factors from
   *  PRD43, 2939 for the decay of spin-1/2 baryons containing bottom and charm
   *  quarks.
   *
   * @see BaryonFormFactor
   * 
   */

class SingletonFormFactor: public BaryonFormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline SingletonFormFactor();

  /**
   * Copy constructor
   */
  inline SingletonFormFactor(const SingletonFormFactor &);

  /**
   * Destructor
   */
  virtual ~SingletonFormFactor();
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
  //@}

  /**
   *  Output the information required to reproduce the object for the particle
   *  properties database
   */
  virtual void dataBaseOutput(ofstream&);
  // output the information for the database

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
  static ClassDescription<SingletonFormFactor> initSingletonFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  SingletonFormFactor & operator=(const SingletonFormFactor &);

private:

  /**
   * The charm quark mass.
   */
  Energy _mcharm;

  /**
   * The strange quark mass.
   */
  Energy _mstrange;

  /**
   *  The mixing angle for the \f$\Lambda\f$.
   */
  double _thetalambda;

  /**
   *  The mixing angle for the \f$\Sigma\f$.
   */
  double _thetasigma;

  /**
   *  The mixing angle for the \f$\Xi\f$.
   */
  double _thetaxi;

  /**
   *  The mixing angle for the \f$\Xi'\f$.
   */
  double _thetaxip;

  /**
   *  The pole masses for the \f$q^2\f$ dependence of the form factors.
   */
  vector<Energy> _polemass;

  /**
   *  The \f$\xi\f$ parameter for the form factors.
   */
  vector<double> _xi;

  /**
   *  The normalisation factor, \f$N_{mM}\f$, for the form factors.
   */
  vector<double> _NmM;

  /**
   *  The mass of the quark for the form factor.
   */
  vector<Energy> _mquark;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * BaryonThreeQuarkModelFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::SingletonFormFactor,1> {
  /** Typedef of the base class of BaryonThreeQuarkModelFormFactor. */
  typedef Herwig::BaryonFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BaryonThreeQuarkModelFormFactor class.
 */
template <>
 struct ClassTraits<Herwig::SingletonFormFactor>
  : public ClassTraitsBase<Herwig::SingletonFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig++::SingletonFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).*/
  static string library() { return "libHwFormFactor.so"; }
};

}

#include "SingletonFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SingletonFormFactor.tcc"
#endif

#endif /* HERWIG_SingletonFormFactor_H */
