// -*- C++ -*-
#ifndef HERWIG_ISGWFormFactor_H
#define HERWIG_ISGWFormFactor_H
//
// This is the declaration of the <ISGWFormFactor class.
//
#include "ScalarFormFactor.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ISGWFormFactor.fh"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The ISGWFormFactor class is the implementation of the ISGW model of 
   *  Phys. Rev. D39, 799 (1989) for the form-factors. It inherits from the 
   *  ScalarFormFactor class and members
   *  for the calculation of the relevant form factors.
   * 
   * @see ScalarFormFactor
   * @see ISGW2FormFactor
   */

class ISGWFormFactor: public ScalarFormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  inline ISGWFormFactor();

  /**
   * Copy constructor
   */
  inline ISGWFormFactor(const ISGWFormFactor &);

  /**
   * Destructor
   */
  virtual ~ISGWFormFactor();
  //@}

public:

  /** @name Form-Factors */
  //@{
  /**
   * The form factor for the weak decay of a scalar to a scalar. 
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f0 The form-factor \f$f_0\f$. 
   * @param fp The form-factor \f$f_+\f$.
   */
  virtual void ScalarScalarFormFactor(Energy2 q2,int iloc,int id0,int id1,Energy m0,
				      Energy m1,Complex & f0,Complex & fp) const;

  /**
   * The form factor for the weak decay of a scalar to a vector.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param V  The form-factor \f$V\f$
   * @param A0 The form-factor \f$A_0\f$
   * @param A1 The form-factor \f$A_1\f$
   * @param A2 The form-factor \f$A_2\f$
   */
  virtual void ScalarVectorFormFactor(Energy2 q2, int iloc, int id0, int id1,
				      Energy m0, Energy m1, Complex & V,
				      Complex & A0,Complex & A1,Complex & A2) const;

  /**
   * The form factor for the weak decay of a scalar to a tensor.
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param h  The form-factor \f$h\f$.
   * @param k  The form-factor \f$k\f$. 
   * @param bp The form-factor \f$b_+\f$.
   * @param bm The form-factor \f$b_-\f$.
   */
  virtual void ScalarTensorFormFactor(Energy2 q2,int iloc,int id0,int id1,Energy m0,
				      Energy m1, Complex & h,Complex & k,
				      Complex & bp, Complex & bm) const;
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

  /** The member which implements all the different form-factors
   * @param q2 The scale \f$q^2\f$.
   * @param iloc The location in the form-factor list.
   * @param id0 The PDG code of the incoming meson.
   * @param id1 The PDG code of the outgoing meson.
   * @param m0 The mass of the incoming meson.
   * @param m1 The mass of the outgoing meson.
   * @param f1 The first  form-factor.
   * @param f2 The second form-factor. 
   * @param f3 The third  form-factor.
   * @param f4 The fourth form-factor.
   */
  void formFactor(Energy2 q2,int iloc,int id0,int id1,Energy m0,
		  Energy m1,Complex & f1,Complex & f2,
		  Complex & f3,Complex & f4) const;
  // general member to calculate all the form-factors

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

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
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ISGWFormFactor> initISGWFormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  ISGWFormFactor & operator=(const ISGWFormFactor &);

private:

  /**
   * The relativistic compension factor
   */
  double _kappa;

  /** @name Quark masses */
  //@{
  /**
   * The down quark mass
   */
  Energy _mdown;

  /**
   * The up quark mass
   */
  Energy _mup;

  /**
   * The strange quark mass
   */
  Energy _mstrange;

  /**
   * The charm quark mass
   */
  Energy _mcharm;

  /**
   * The bottom quark mass
   */
  Energy _mbottom;

  /**
   * The masses of the quarks as a vector
   */
  vector<Energy> _mquark;
  //@}


  /** @name Wave function parameters */
  //@{
  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{d}\f$ 
   */
  Energy _betaSud;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{s}\f$ 
   */
  Energy _betaSus;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{c}\f$ 
   */
  Energy _betaSuc;

  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{b}\f$ 
   */
  Energy _betaSub;

  /**
   * The s-wave variational parameters as a vector.
   */
  vector<Energy> _betaS;

  /**
   * The wavefunction p-wave \f$\beta\f$ variational parameters for  \f$u\bar{d}\f$ 
   */
  Energy _betaPud;
  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{s}\f$ 
   */
  Energy _betaPus;
  /**
   * The wavefunction s-wave \f$\beta\f$ variational parameters for  \f$u\bar{c}\f$ 
   */
  Energy _betaPuc;

  /**
   * The p-wave variational parameters as a vector
   */
  vector<Energy> _betaP;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/**
 * This template specialization informs ThePEG about the base class of
 * ISGWFormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::ISGWFormFactor,1> {
  /** Typedef of the base class of ISGWFormFactor. */
  typedef Herwig::ScalarFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ISGWFormFactor class.
 */
template <>
struct ClassTraits<Herwig::ISGWFormFactor>
  : public ClassTraitsBase<Herwig::ISGWFormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig++::ISGWFormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwFormFactors.so"; }
};

}

#include "ISGWFormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ISGWFormFactor.tcc"
#endif

#endif /* HERWIG_ISGWFormFactor_H */
