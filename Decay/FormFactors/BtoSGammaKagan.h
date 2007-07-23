// -*- C++ -*-
#ifndef HERWIG_BtoSGammaKagan_H
#define HERWIG_BtoSGammaKagan_H
//
// This is the declaration of the BtoSGammaKagan class.
//

#include "Herwig++/Utilities/Interpolator.h"
#include "BtoSGammaHadronicMass.h"
#include "ThePEG/Config/Complex.h"
#include "ThePEG/Config/Constants.h"
#include "BtoSGammaKagan.fh"

namespace Herwig {

using namespace ThePEG;
using namespace ThePEG::Constants;

/** \ingroup Decay
 *
 * The BtoSGammaKagan class implements the model of hep-ph/9805303 for the 
 * hadronic mass spectrum in \f$b\to s \gamma\f$ decays.
 */
class BtoSGammaKagan: public BtoSGammaHadronicMass {

  /**
   * Class for the integration is a friend to access private members
   */
  friend class KaganIntegrand;

public:

  /**
   * The default constructor.
   */
  inline BtoSGammaKagan();

  /**
   * Returns the hadronic mass.
   * @param mb The mass of the decaying B meson
   * @param mquark The minimum mass of the hadronic system based on the consistuent quark
   * masses.
   * @return The hadronic mass
   */
  virtual Energy hadronicMass(Energy mb,Energy mquark);

  /**
   * Output the setup information for the particle database
   * @param os The stream to output the information to
   * @param header Whether or not to output the information for MySQL
   * @param create Whether or not to add a statement creating the object
   */
  virtual void dataBaseOutput(ofstream & os,bool header,bool create) const;

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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:

  /**
   *  Members which return integrands
   */
  //@{
  /**
   *  Operator to return the integrand for the \f$s_{22}(y)\f$ function
   *  or \f$s_{27}(y)\f$ functions of hep-ph/9805303
   *  depending on the value of _iopt to be integrated
   */
  inline double operator ()(double argument) const;
  typedef double ValType;
  typedef double ArgType;

  /**
   *  Operator to return the integrand of the smeared function or
   *  Fermi function depending on the value of _iopt to be integrated
   */
  inline InvEnergy operator ()(Energy argument) const;
  //@}

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
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  inline virtual void doinitrun();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<BtoSGammaKagan> initBtoSGammaKagan;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  BtoSGammaKagan & operator=(const BtoSGammaKagan &);

private:

  /** @name Functions to calculate the mass spectrum */
  //@{
  /**
   * The derivative of the Sudakov form-factor from hep-ph/9805303
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   * @param alphaS The strong coupling, \f$\alpha_S\f$.
   */
  inline double Delta(double y, double alphaS) const;

  /**
   * Kinematic function from semi-leptonic decay for normaalisation
   */
  inline double semiLeptonicf() const;

  /**
   *  \f$s_{22}(y)\f$ function from hep-ph/9805303. Due to the integration
   * required this function is computed by interpolation.
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   */
  inline double s22(double y) const;

  /**
   *  \f$s_{27}(y)\f$ function from hep-ph/9805303. Due to the integration
   * required this function is computed by interpolation.
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   */
  inline double s27(double y) const;

  /**
   *  \f$s_{77}(y)\f$ function from hep-ph/9805303
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   */
  inline double s77(double y) const;

  /**
   *  \f$s_{78}(y)\f$ function from hep-ph/9805303
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   */
  inline double s78(double y) const;

  /**
   *  \f$s_{88}(y)\f$ function from hep-ph/9805303
   * @param y Ratio \f$E_\gamma/E^{\rm max}_\gamma\f$.
   */
  inline double s88(double y) const;

  /**
   *  The real part of the \f$G(t)\f$ function from hep-ph/9805303
   */
  inline double realG(double) const;

  /**
   *  The imaginary part of the \f$G(t)\f$ function from hep-ph/9805303
   */
  inline double imagG(double) const;

  /**
   *  Strong coupling \f$\alpha_S\f$ at the scale \f$Q\f$
   * @param Q The scale.
   */
  inline double alphaS(Energy Q);

  /**
   *   Calculate the wilson coefficients we need
   */
  void calculateWilsonCoefficients();

  /**
   * The \f$K'_{NLO}(1-y)\f$ function at parton level from hep-ph/9805303
   */
  inline double KNLO(double) const;
  //@}

private:


  /**
   *  Initialisation of mass spectrum
   */
  bool _initialize;

  /**
   *  Quark masses and related parameters
   */
  //@{
  /**
   *  The top quark mass
   */
  Energy _mt;

  /**
   * bottom quark mass
   */
  Energy _mb;

  /**
   * charm quark mass
   */
  Energy _mc;

  /**
   * strange quark mass
   */
  Energy _ms;

  /**
   *  Ratio of the strange quark mass to the bottom quark mass
   */
  double _msovermb;

  /**
   * The ratio of the charm to bottom quark masses squared, \f$(m_c/m_b)^2\f$
   */
  double _zratio;
  //@}

  /**
   * The hadronic \f$\lambda_2\f$ parameter from hep-ph/9805303.
   */
  Energy2 _lambda2;

  /**
   *  Masses of other particles
   */
  //@{
  /**
   *  The W mass
   */
  Energy _mw;

  /**
   *  the Z mass
   */
  Energy _mz;

  /**
   *  Mass of the decaying B meson.
   */
  Energy _MB;
  //@}

  /** @name Wilson coefficients, couplings and \f$\beta\f$-function coefficients*/
  //@{
  /**
   * The leading order \f$c_2\f$ coefficient.
   */
  double _c20;

  /**
   * The leading order \f$c_7\f$ coefficient.
   */
  double _c70;

  /**
   * The leading order \f$c_8\f$ coefficient.
   */
  double _c80;

  /**
   *  First \f$\beta\f$-function coefficient
   */
  double _beta0;

  /**
   *  Second \f$\beta\f$-function coefficient
   */
  double _beta1;

  /**
   *  The electromagentic coupling
   */
  double _alpha;

  /**
   *  The strong coupling at the Z mass
   */
  double _alphaSZ;

  /**
   *  The renormalisation scale
   */
  Energy _mub;

  /**
   *  the strong coupling at the renormalisation scale \f$\mu_b\f$.
   */
  double _alphaSM;

  /**
   *   The CKM perfactor for the decay
   */
  double _ckm;
  /**
   *  Pre-factor for the correction term involving \f$\Delta(y)\f$.
   */
  double _delta;
  //@}

  /**
   *  Interpolators for the integrate functions and related parameters
   */
  //@{
  /**
   *  Interpolator for the \f$s_{22}\f$ function
   */
  Interpolator<double,double>::Ptr _s22inter;

  /**
   *  Interpolator for the \f$s_{27}\f$ function
   */
  Interpolator<double,double>::Ptr _s27inter;

  /**
   *  Interpolator for the spectrum
   */
  Interpolator<InvEnergy,Energy>::Ptr _pmHinter;

  /**
   *  Values of \f$m_H\f$ for the interpolation of the spectrum
   */
  vector<Energy> _mHinter;

  /**
   *  Values of the differential rate for the interpolation of the spectrum
   */
  vector<InvEnergy> _spectrum;

  /**
   *  Maximum value of the spectrum for unweighting
   */
  InvEnergy _spectmax;

  /**
   *  Maximum number of tries for unweighting
   */
  unsigned int _maxtry;
  //@}

  /**
   *  Parameters for the Fermi function
   */
  //@{
  /**
   *  The \f$\bar{\Lambda}\f$ parameter from hep-ph/9805303.
   */
  Energy _fermilambda;

  /**
   *  The power from  from hep-ph/9805303.
   */
  double _fermia;

  /**
   * The normalisation from hep-ph/9805303.
   */
  InvEnergy _ferminorm;

  /**
   * \f$\lambda_1\f$ scale related to the kinetic energy of the b quark.
   */
  Energy2 _fermilambda1;
  //@}

  /**
   *   Techincal parameters for the integration of the spectrum
   */
  //@{
  /**
   * Cut-off parameter to avoid the singularity at y=1
   */
  double _ycut;

  /**
   * Value of the energy fraction for which the integral is being performed
   */
  double _y;

  /**
   *  Cut-off on the photon energies
   */
  double _deltacut;

  /**
   * Number of points for the interpolation of the s functions
   */
  unsigned int _nsfunct;

  /**
   * Number of points for the interpolation of the spectrum
   */
  unsigned int _nspect;

  /**
   *  The function currently being integrated
   */
  unsigned int _iopt;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of BtoSGammaKagan. */
template <>
struct BaseClassTrait<Herwig::BtoSGammaKagan,1> {
  /** Typedef of the first base class of BtoSGammaKagan. */
  typedef Herwig::BtoSGammaHadronicMass NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the BtoSGammaKagan class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::BtoSGammaKagan>
  : public ClassTraitsBase<Herwig::BtoSGammaKagan> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::BtoSGammaKagan"; }
  /** Return the name of the shared library be loaded to get
   *  access to the BtoSGammaKagan class and every other class it uses
   *  (except the base class). */
  static string library() { return "HwFormFactors.so"; }
};

/** @endcond */

}

namespace Herwig {

/**
 *  A struct for the integrand which can access the dimensional value
 *  member of the BtoSGammaKagan class
 */
struct KaganIntegrand {

  /**
   *  The constructor
   */
  inline KaganIntegrand(BtoSGammaKaganPtr in) : _kagan(in) {};

  /**
   * Get the function value
   */
  inline InvEnergy operator ()(Energy arg) const {return (*_kagan)(arg);}
  typedef InvEnergy ValType;
  typedef Energy    ArgType;

  /**
   *  A pointer to the form factor to supply the integrand.
   */
  BtoSGammaKaganPtr _kagan;
};
}

#include "BtoSGammaKagan.icc"

#endif /* HERWIG_BtoSGammaKagan_H */
