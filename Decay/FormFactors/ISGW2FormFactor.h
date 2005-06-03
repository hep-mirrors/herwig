// -*- C++ -*-
#ifndef HERWIG_ISGW2FormFactor_H
#define HERWIG_ISGW2FormFactor_H
//
// This is the declaration of the ISGW2FormFactor class.
//
#include "ScalarFormFactor.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ISGW2FormFactor.fh"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {
using namespace ThePEG;

  /** \ingroup Decay
   *
   *  The ISGW2FormFactor class is the implementation of 
   *  the ISGW2 model of Phys. Rev. D52, 2783 (1995) for the scalar meson form
   *  factors.
   *
   *  It inherits from the ScalarFormFactor class and implements
   *  the calculation of the relevant form factors.
   *
   * @see ScalarFormFactor
   * @see ISGWFormFactor
   */

class ISGW2FormFactor: public ScalarFormFactor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  ISGW2FormFactor();

  /**
   * Copy constructor
   */
  inline ISGW2FormFactor(const ISGW2FormFactor &);

  /**
   * Destructor
   */
  virtual ~ISGW2FormFactor();
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
  virtual void ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1,Complex & f0,
				      Complex & fp) const;

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
  virtual void ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1,
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
  virtual void ScalarTensorFormFactor(Energy2 q2,unsigned int iloc,int id0,int id1,
				      Energy m0,Energy m1, Complex & h,Complex & k,
				      Complex & bp, Complex & bm) const;
  //@}

  /**
   * Output the setup information for the particle database
   */
  virtual void dataBaseOutput(ofstream &);

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
  void formFactor(Energy2 q2,unsigned int iloc,int id0,int id1,Energy m0,
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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
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
  static ClassDescription<ISGW2FormFactor> initISGW2FormFactor;

  /**
   * Private and non-existent assignment operator.
   */
  ISGW2FormFactor & operator=(const ISGW2FormFactor &);

protected:

  /**
   * The saturated \f$\alpha_S\f$ used to calculate the form-factors.
   * @param mass Mass scale to work out the number of flavours.
   * @param q2 \f$q^2\f$ the scale.
   * @return the value of \f$\alpha_S\f$.
   */
  double alphaS(Energy mass, Energy2 q2) const;

private:

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
  Energy _mquark[5];
  //@}

  /** @name Wave function parameters for the \f$1^1S_0\f$ level.*/
  //@{

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$u\bar{d}\f$ 
   */
  Energy _beta1S0ud;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$u\bar{s}\f$ 
   */
  Energy _beta1S0us;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$s\bar{s}\f$ 
   */
  Energy _beta1S0ss;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$c\bar{u}\f$ 
   */
  Energy _beta1S0cu;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$c\bar{s}\f$ 
   */
  Energy _beta1S0cs;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$u\bar{b}\f$ 
   */
  Energy _beta1S0ub;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$s\bar{b}\f$ 
   */
  Energy _beta1S0sb;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$c\bar{c}\f$ 
   */
  Energy _beta1S0cc;

  /**
   * The wavefunction \f$1^1S_0\f$ \f$\beta\f$ variational parameter for \f$b\bar{c}\f$ 
   */
  Energy _beta1S0bc;

  /**
   *  The wavefunction parameters as an array
   */
  Energy _beta1S0[5][5];

  /**
   *  The masses as a array
   */
  Energy _mass1S0[5][5];
  //@}

  /** @name Wave function parameters for the \f$1^3S_1\f$ level.*/
  //@{
  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$u\bar{d}\f$ 
   */
  Energy _beta3S1ud;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$u\bar{s}\f$ 
   */
  Energy _beta3S1us;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$s\bar{s}\f$ 
   */
  Energy _beta3S1ss;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$c\bar{u}\f$ 
   */
  Energy _beta3S1cu;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$c\bar{s}\f$ 
   */
  Energy _beta3S1cs;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$u\bar{b}\f$ 
   */
  Energy _beta3S1ub;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$s\bar{b}\f$ 
   */
  Energy _beta3S1sb;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$c\bar{c}\f$ 
   */
  Energy _beta3S1cc;

  /**
   * The wavefunction \f$1^3S_1\f$ \f$\beta\f$ variational parameter for \f$b\bar{c}\f$ 
   */
  Energy _beta3S1bc;

  /**
   * The wavefunction paramaeters as an array.
   */
  Energy _beta3S1[5][5];
  //@}


  /** @name Wave function parameters for the \f$1P\f$ levels.*/
  //@{

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$u\bar{d}\f$ 
   */
  Energy _beta1Pud;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$u\bar{s}\f$ 
   */
  Energy _beta1Pus;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$s\bar{s}\f$ 
   */
  Energy _beta1Pss;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$c\bar{u}\f$ 
   */
  Energy _beta1Pcu;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$c\bar{s}\f$ 
   */
  Energy _beta1Pcs;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$u\bar{b}\f$ 
   */
  Energy _beta1Pub;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$s\bar{b}\f$ 
   */
  Energy _beta1Psb;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$c\bar{c}\f$ 
   */
  Energy _beta1Pcc;

  /**
   * The wavefunction \f$1P\f$ \f$\beta\f$ variational parameter for \f$b\bar{c}\f$ 
   */
  Energy _beta1Pbc;

  /**
   * The wavefunction paramaeters as an array.
   */
  Energy _beta1P[5][5];

  /**
   * The spin-1/2 masses
   */
  // the 1/2 spin masses
  Energy _massPoh[5][5];

  /**
   * The spin-3/2 masses
   */
  Energy _massPth[5][5];
  //@}

  /**@name Parameters for the strong coupling*/
  //@{
  /**
   *  The cut-off value of \f$\alpha_S\f$.
   */
  double _alphamuQM;
  /**
   * The values of \f$\alpha_S\f$ at the quark masses.
   */
  double _alphaQ[5];
  //@}

  /**@name Relativistic correction factors */
  //@{
  /**
   * The correction factor for \f$D\to\rho\f$.
   */
  double _CfDrho;

  /**
   * The correction factor for \f$D\to K^*\f$.
   */
  double _CfDKstar;

  /**
   * The correction factor for \f$D_s\to K^*\f$.
   */
  double _CfDsKstar;

  /**
   * The correction factor for \f$D_s\to\phi\f$.
   */
  double _CfDsphi;

  /**
   * The correction factor for \f$B\to\rho\f$.
   */
  double _CfBrho;

  /**
   * The correction factor for \f$B\to D^*\f$.
   */
  double _CfBDstar;

  /**
   * The correction factor for \f$B_s\to K^*\f$.
   */
  double  _CfBsKstar;

  /**
   * The correction factor for \f$B_s\to D^*\f$.
   */
  double _CfBsDstar;

  /**
   * The correction factor for \f$B_c\to D^*\f$.
   */
  double _CfBcDstar;

  /**
   * The correction factor for \f$B_c\to\psi\f$.
   */
  double _CfBcpsi;

  /**
   * The correction factor for \f$B_c\to B_s^*\f$.
   */
  double _CfBcBsstar;

  /**
   * The correction factor for \f$B_c\to B^*\f$.
   */
  double _CfBcBstar;
  //@}

  /**
   * The \f$\eta-\eta'\f$ mixing angle 
   */
  double _thetaeta;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {
/**
 * This template specialization informs ThePEG about the base class of
 * ISGW2FormFactor.
 */
template <>
 struct BaseClassTrait<Herwig::ISGW2FormFactor,1> {
  /** Typedef of the base class of ISGW2FormFactor. */
   typedef Herwig::ScalarFormFactor NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ISGW2FormFactor class.
 */
template <>
 struct ClassTraits<Herwig::ISGW2FormFactor>
  : public ClassTraitsBase<Herwig::ISGW2FormFactor> {
  /** Return the class name. */
  static string className() { return "Herwig++::ISGW2FormFactor"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class and every other class it uses
   * (except the base class).
   */
  static string library() { return "libHwFormFactor.so"; }
};

}

#include "ISGW2FormFactor.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ISGW2FormFactor.tcc"
#endif

#endif /* HERWIG_ISGW2FormFactor_H */

