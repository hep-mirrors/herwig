// -*- C++ -*-
#ifndef HERWIG_TopDecayMECorrection_H
#define HERWIG_TopDecayMECorrection_H
//
// This is the declaration of the TopDecayMECorrection class.
//

#include "MECorrectionBase.h"
#include "TopDecayMECorrection.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * The TopDecayMECorrection class implements the matrix element correction
 * for top decay.
 *
 * @see \ref TopDecayMECorrectionInterfaces "The interfaces"
 * defined for TopDecayMECorrection.
 */
class TopDecayMECorrection: public MECorrectionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline TopDecayMECorrection();

  /**
   * The copy constructor.
   */
  inline TopDecayMECorrection(const TopDecayMECorrection &);

  /**
   * The destructor.
   */
  virtual ~TopDecayMECorrection();
  //@}

public:

  /**
   *  Members to override those in the base class and implemented 
   *  the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   */
  virtual bool canHandle(ShowerTreePtr);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,Branching br);
  //@}

private:

  /**
   *  Apply the hard matrix element
   */
  vector<Lorentz5Momentum> applyHard(const ParticleVector &p,double,double);

  /**
   *  Get the weight for hard emission
   */
  double getHard(double, double);

  /**
   *  This function is auxiliary to the function $x_{a}$ (hXAB).
   */
  inline double XGBR1(double,double,double,int);
  /**
   *  This function is auxiliary to the function $x_{a}$ (hXAB).
   */
  inline double XGBR2(double,double,double,int);
  /**
   *  This function is auxiliary to the function $x_{a}$ (hXAB).
   */
  inline double KTR(double,double,double,double,int);
  /**
   *  This function determines $x_{a}$ as a function of $x_{g}$ 
   *  and $\kappa$ where $\kappa$ pertains to emissions from the 
   *  b.
   */
  inline double XAB(double,double,double,double,double,int);
  /**
   *  This function determines the point ($x_{g}$) where the condition that 
   *  $x_{a}$ be real supersedes that due to the external input 
   *  $\tilde{\kappa}$ where, again, $\kappa$ pertains to emissions from the 
   *  b.
   */
  inline double XGBCUT(double,double,double);
  /**
   *  This function determines the minimum value of $x_{a}$ 
   *  for a given $\tilde{\kappa}$ where $\kappa$ pertains to
   *  emissions from the c.
   */
  inline double XACCUT(double,double,double);
  /**
   *  This function is auxiliary to the function $x_{g}$ (hXGC).
   */
  inline double Z(double,double,double,double,double,int,int); 
  /**
   *  This function determines $x_{g}$ as a function of $x_{a}$ 
   *  and $\kappa$ where $\kappa$ pertains to emissions from the 
   *  c. It is multivalued, one selects a branch according to the
   *  second to last integer flag (+/-1). The last integer flag
   *  is used to select whether (1) or not (0) you wish to have the 
   *  function for the special case of the full phase space, in which
   *  case the fifth argument $\kappa$ is irrelevant.
   */
  inline double XGC(double,double,double,double,double,int,int); 
  /**
   *  This function, $x_{g,c=0}^{-1}$, returns $x_{a}$ as a function 
   *  of $x_{g}$ for the special case of c=0, for emissions from c 
   *  (the b-quark). The third input is $\tilde{\kappa}$ which pertains 
   *  to emissions from c.
   */
  inline double XGINVC0(double,double,double,double); 
  /**
   *  For a given value of $x_{g}$ this returns the maximum value of $x_{a}$  
   *  in the dead region.
   */
  inline double APPROXDEADMAXXA(double,double,double); 
  /**
   *  For a given value of $x_{g}$ this returns the maximum value of $x_{a}$  
   *  in the dead region.
   */
  inline double APPROXDEADMINXA(double,double,double); 
  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are in the allowed region, the kinematically accessible phase 
   *  space.
   */
  inline bool INTHEALLOWEDREGION(double,double,double,double,double); 
  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the approximate dead region.
   */
  inline bool INTHEAPPROXDEADREGION(double,double,double,double,
                                    double,double,double); 
  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the dead region.
   */
  inline bool INTHEDEADREGION(double,double,double,double,
                              double,double,double); 

  /**
   *  This function returns values of ($x_{g}$,$x_{a}$) distributed 
   *  according to $\left(1+a-x_{a}\right)^{-1}x_{g}^{-2}$ in the 
   *  approximate dead region.  
   */
  inline double DEADREGIONXGXA(double,double); 
  /**
   *  This is similar to the subroutine of the same name (stolen) in
   *  herwig6507. It takes a 5-momentum and returns a rotation matrix 
   *  such that which acts on the input 5-momentum so as to rotate
   *  it to point in the +Z direction.  
   */
  inline HepLorentzRotation HWUROT(Lorentz5Momentum); 
  /**
   *  This returns a random rotation about the +Z direction.  
   */
  inline HepLorentzRotation RANDOMZROTATION(); 
  /**
   *  This is the same as HWUSQR in herwig6507. It returns the square 
   *  root of the absolute value of the input, multiplied by the sign
   *  of the input.  
   */
  inline double HWUSQR(double); 
  /**
   *  This routine can be removed, it is being used in dofinish
   *  as part of debugging, to print topdrawer files.
   */ 
  inline int KPRINTER(double,double);  
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

  /**
   *  Full matrix element with a factor of \f$\frac{\alpha_SC_F}{x_g^2\pi}\f$ removed.
   * @param xw The momentum fraction of the W boson
   * @param xg The momentum fraction of the gluon.
   */
  double me(double xw, double xg);

  /**
   * Full matrix element with a factor of 
   * \f$\frac{\alpha_{S}C_{F}}{\left(1.+A-x_{W}\right)x_{g}^{2}\pi}\f$ removed.
   * This is the expression inside the braces in equation (6.31) of
   * Gieseke,Stephens & Webber, JHEP12(2003) 045.
   * @param A  The squared ratio of the b-quark mass to the top quark mass.
   * @param XA The momentum fraction of the W boson.
   * @param XG The momentum fraction of the gluon.
   */
  inline double BRACES(double XG, double XA);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<TopDecayMECorrection> initTopDecayMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  TopDecayMECorrection & operator=(const TopDecayMECorrection &);

private:

  /**
   *  The mass of the W boson
   */
  Energy _ma;

  /**
   *  The mass of the bottom quark
   */
  Energy _mc;

  /**
   *  The top mass
   */
  Energy _mt;

  /**
   *  The gluon mass.
   */
  Energy _mg;

  /**
   *  The mass ratio for the W.
   */
  double _a;

  /**
   *  The mass ratio for the bottom.
   */
  double _c;

  /**
   *  The mass ratio for the gluon.
   */
  double _g;

  /**
   *  Two times the energy fraction of a.
   */
  double _ktb;

  /**
   *  Two times the energy fraction of the gluon.
   */
  double _ktc;

  /**
   *  Two times the energy fraction of the gluon.
   */
  double _xg;

  /**
   *  Two times the energy fraction of a.
   */
  double _xa;

  /**
   *  Two times the energy fraction of c.
   */
  double _xc;

  /**
   *  This determines the hard matrix element importance 
   *  sampling in _xg. _xg_sampling=2.0 samples as 1/xg^2.
   */
  double _xg_sampling;

  /**
   *  The enhancement factor for initial-state radiation
   */
  double _initialenhance;

  /**
   *  The enhancement factor for final-state radiation
   */
  double _finalenhance;

  /**
   *  This flag determines whether the T2 region in the decay shower
   *  (JHEP12(2003)_045) is populated by the ME correction (true) or
   *  the shower from the decaying particle.
   */
  bool _use_me_for_t2;


};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of TopDecayMECorrection. */
template <>
struct BaseClassTrait<Herwig::TopDecayMECorrection,1> {
  /** Typedef of the first base class of TopDecayMECorrection. */
  typedef Herwig::MECorrectionBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TopDecayMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TopDecayMECorrection>
  : public ClassTraitsBase<Herwig::TopDecayMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig++::TopDecayMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TopDecayMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class TopDecayMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#include "TopDecayMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDecayMECorrection.tcc"
#endif

#endif /* HERWIG_TopDecayMECorrection_H */
