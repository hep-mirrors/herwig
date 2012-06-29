// -*- C++ -*-
#ifndef HERWIG_BaryonSimpleFormFactor_H
#define HERWIG_BaryonSimpleFormFactor_H
//
// This is the declaration of the BaryonSimpleFormFactor class.
//
#include "BaryonFormFactor.h"

namespace Herwig {
using namespace ThePEG;

/** \ingroup Decay
 *
 *  The BaryonSimpleFormFactor class is a simple model for the form-factors
 *  for the semi-leptonic decay of the light (i.e. uds) baryons. The form-factors
 *  are assumed to be constant and are taken from the quark model results
 *  of PRD25, 206 (1982).
 *
 * @ see BaryonFormFactor
 */
class BaryonSimpleFormFactor: public BaryonFormFactor {

public:

  /**
   * Default constructor
   */
  BaryonSimpleFormFactor();

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
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
  static ClassDescription<BaryonSimpleFormFactor> initBaryonSimpleFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  BaryonSimpleFormFactor & operator=(const BaryonSimpleFormFactor &);

private:

  /**
   * The axial vector coupling in \f$\beta\f$ decay.
   */
  double _gA;

  /**
   *  The \f$\frac{D}{D+F}\f$ ratio.
   */
  double _alphaD;

  /**
   *  The \f$\eta_V\f$ \f$SU(3)\f$ breaking parameter.
   */
  double _etaV;

  /**
   *  The \f$\eta_V\f$ \f$SU(3)\f$ breaking parameter.
   */
  double _etaA;

  /**
   *  The \f$\rho_E\f$  \f$SU(3)\f$ breaking parameter for the electic dipole moment.
   */
  double _rhoE;

  /**
   *  The \f$\rho_M\f$  \f$SU(3)\f$ breaking parameter for the electic dipole moment.
   */
  double _rhoM;

  /**
   *  The calculated constants for the \f$f_1\f$ form factor.
   */
  vector<double> _f1;

  /**
   *  The calculated constants for the \f$f_2\f$ form factor.
   */
  vector<double> _f2;

  /**
   *  The calculated constants for the \f$g_1\f$ form factor.
   */
  vector<double> _g1;

  /**
   *  The calculated constants for the \f$g_2\f$ form factor.
   */
  vector<double> _g2;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * BaryonSimpleFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::BaryonSimpleFormFactor,1> {
  /** Typedef of the base class of BaryonSimpleFormFactor. */
  typedef Herwig::BaryonFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BaryonSimpleFormFactor class.
 */
template <>
struct ClassTraits<Herwig::BaryonSimpleFormFactor>
  : public ClassTraitsBase<Herwig::BaryonSimpleFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig::BaryonSimpleFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

#endif /* HERWIG_BaryonSimpleFormFactor_H */
