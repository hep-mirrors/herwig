// -*- C++ -*-
#ifndef HERWIG_LambdabExcitedLambdacSumRuleFormFactor_H
#define HERWIG_LambdabExcitedLambdacSumRuleFormFactor_H
// This is the declaration of the LambdabExcitedLambdacSumRuleFormFactor class.

#include "BaryonFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>LambdabExcitedLambdacSumRuleFormFactor</code> class implements the
 *  form-factors of hep-ph/0012114 for the decay of the \f$Lambda_b^0\f$
 *  to excited \f$\Lambda^+_c\f$ states.
 *
 * @see BaryonFormFactor.
 * 
 */
class LambdabExcitedLambdacSumRuleFormFactor: public BaryonFormFactor {

public:

  /**
   * Default constructor
   */
  LambdabExcitedLambdacSumRuleFormFactor();

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LambdabExcitedLambdacSumRuleFormFactor> initLambdabExcitedLambdacSumRuleFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  LambdabExcitedLambdacSumRuleFormFactor & operator=(const LambdabExcitedLambdacSumRuleFormFactor &);

private:

  /**
   * The intercept for the form factor. 
   */
  double _xi1;

  /**
   *  The slope for the form factor.
   */
  double _rho2;
};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of LambdabExcitedLambdacSumRuleFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::LambdabExcitedLambdacSumRuleFormFactor,1> {
  /** Typedef of the base class of  LambdabExcitedLambdacSumRuleFormFactor.*/
  typedef Herwig::BaryonFormFactor NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
struct ClassTraits<Herwig::LambdabExcitedLambdacSumRuleFormFactor>
  : public ClassTraitsBase<Herwig::LambdabExcitedLambdacSumRuleFormFactor> {
  /** Return the class name. */
  static string className() {return "Herwig::LambdabExcitedLambdacSumRuleFormFactor";}
  /**
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwFormFactors.so"; }

};

/** @endcond */

}

#endif /* HERWIG_LambdabExcitedLambdacSumRuleFormFactor_H */
