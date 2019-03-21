// -*- C++ -*-
#ifndef HERWIG_LightBaryonQuarkModelFormFactor_H
#define HERWIG_LightBaryonQuarkModelFormFactor_H
// This is the declaration of the LightBaryonQuarkModelFormFactor class.

#include "BaryonFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The <code>LightBaryonQuarkModelFormFactor</code> class implements the 
 *  quark model calculation of hep-ph/9409272 for the form-factors for 
 *  the decay of baryons containing the light quarks. This can be adapted for
 *  other models.
 *
 * @see BaryonFormFactor
 * 
 */
class LightBaryonQuarkModelFormFactor: public BaryonFormFactor {

public:

  /**
   * Default constructor
   */
  LightBaryonQuarkModelFormFactor();

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
  //@}

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving and
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LightBaryonQuarkModelFormFactor> initLightBaryonQuarkModelFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  LightBaryonQuarkModelFormFactor & operator=(const LightBaryonQuarkModelFormFactor &) = delete;

private:

  /**
   * The form-factor \f$f_1\f$ at \f$q^2=0\f$.
   */
  vector<double> _f1;

  /**
   * The form-factor \f$f_2\f$ at \f$q^2=0\f$.
   */
  vector<InvEnergy> _f2;

  /**
   * The form-factor \f$g_1\f$ at \f$q^2=0\f$.
   */
  vector<double> _g1;

  /**
   * The form-factor \f$g_2\f$ at \f$q^2=0\f$.
   */
  vector<InvEnergy> _g2;

  /**
   * The mass for the energy dependence of the \f$f_1\f$ form factor.
   */
  vector<Energy> _Lambdaf1;

  /**
   * The mass for the energy dependence of the \f$f_2\f$ form factor.
   */
  vector<Energy> _Lambdaf2;

  /**
   * The mass for the energy dependence of the \f$g_1\f$ form factor.
   */
  vector<Energy> _Lambdag1;

  /**
   * The mass for the energy dependence of the \f$g_2\f$ form factor.
   */
  vector<Energy> _Lambdag2;

};

}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of LightBaryonQuarkModelFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::LightBaryonQuarkModelFormFactor,1> {
  /** Typedef of the base class of LightBaryonQuarkModelFormFactor. */
   typedef Herwig::BaryonFormFactor NthBase;
};

template <>
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
 struct ClassTraits<Herwig::LightBaryonQuarkModelFormFactor>
  : public ClassTraitsBase<Herwig::LightBaryonQuarkModelFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::LightBaryonQuarkModelFormFactor"; }
  /**
   * Return the class name.
   * Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwFormFactors.so"; }

};

/** @endcond */

}

#endif /* HERWIG_LightBaryonQuarkModelFormFactor_H */
