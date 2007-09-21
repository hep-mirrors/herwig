// -*- C++ -*-
#ifndef HERWIG_TopDecayMECorrection_H
#define HERWIG_TopDecayMECorrection_H
//
// This is the declaration of the TopDecayMECorrection class.
//

#include "QTildeMECorrection.h"
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
class TopDecayMECorrection: public QTildeMECorrection {

public:

  /**
   * The default constructor.
   */
  inline TopDecayMECorrection();

  /**
   *  Members to override those in the base class and implemented 
   *  the matrix element correction
   */
  //@{
  /**
   *  Can the matrix element correction handle a given hard process or decay
   * @param tree The shower tree currently being showered
   * @param initial The initial-state radiation enhancement factor
   * @param final   The final-state radiation enhancement factor
   * @param evolver Pointer to the Evolver.
   */
  virtual bool canHandle(ShowerTreePtr tree,double & initial,
			 double & final, EvolverPtr evolver);

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual void applyHardMatrixElementCorrection(ShowerTreePtr);

  /**
   * Apply the soft matrix element correction
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
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
   *  This function is auxiliary to the function \f$x_{a}\f$ (hXAB).
   */
  inline double xgbr(int);

  /**
   *  This function is auxiliary to the function \f$x_{a}\f$ (hXAB).
   */
  inline double ktr(double,int);

  /**
   *  This function determines \f$x_{a}\f$ as a function of \f$x_{g}\f$ 
   *  and \f$\kappa\f$ where \f$\kappa\f$ pertains to emissions from the 
   *  b.
   */
  inline double xab(double,double,int);

  /**
   *  This function determines the point (\f$x_{g}\f$) where the condition that 
   *  \f$x_{a}\f$ be real supersedes that due to the external input 
   *  \f$\tilde{\kappa}\f$ where, again, \f$\kappa\f$ pertains to emissions from the 
   *  b.
   */
  inline double xgbcut(double);

  /**
   *  This function determines the minimum value of \f$x_{a}\f$ 
   *  for a given \f$\tilde{\kappa}\f$ where \f$\kappa\f$ pertains to
   *  emissions from the c.
   */
  inline double xaccut(double);

  /**
   *  This function is auxiliary to the function \f$x_{g}\f$ (hXGC).
   */
  inline double z(double,double,int,int); 

  /**
   *  This function determines \f$x_{g}\f$ as a function of \f$x_{a}\f$ 
   *  and \f$\kappa\f$ where \f$\kappa\f$ pertains to emissions from the 
   *  c. It is multivalued, one selects a branch according to the
   *  second to last integer flag (+/-1). The last integer flag
   *  is used to select whether (1) or not (0) you wish to have the 
   *  function for the special case of the full phase space, in which
   *  case the fifth argument \f$\kappa\f$ is irrelevant.
   */
  inline double xgc(double,double,int,int); 

  /**
   *  This function, \f$x_{g,c=0}^{-1}\f$, returns \f$x_{a}\f$ as a function 
   *  of \f$x_{g}\f$ for the special case of c=0, for emissions from c 
   *  (the b-quark). The third input is \f$\tilde{\kappa}\f$ which pertains 
   *  to emissions from c.
   */
  inline double xginvc0(double,double); 

  /**
   *  For a given value of \f$x_{g}\f$ this returns the maximum value of \f$x_{a}\f$  
   *  in the dead region.
   */
  inline double approxDeadMaxxa(double,double,double); 

  /**
   *  For a given value of \f$x_{g}\f$ this returns the maximum value of \f$x_{a}\f$  
   *  in the dead region.
   */
  inline double approxDeadMinxa(double,double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are in the allowed region, the kinematically accessible phase 
   *  space.
   */
  inline bool inTheAllowedRegion(double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the approximate dead region.
   */
  inline bool inTheApproxDeadRegion(double,double,
                                    double,double); 

  /**
   *  This function returns true or false according to whether the values
   *  xg,xa are exactly in the dead region.
   */
  inline bool inTheDeadRegion(double,double,
                              double,double); 

  /**
   *  This function returns values of (\f$x_{g}\f$,\f$x_{a}\f$) distributed 
   *  according to \f$\left(1+a-x_{a}\right)^{-1}x_{g}^{-2}\f$ in the 
   *  approximate dead region.  
   */
  inline double deadRegionxgxa(double,double); 

  /**
   *  This rotation takes a 5-momentum and returns a rotation matrix 
   *  such that it acts on the input 5-momentum so as to
   *  make it point in the +Z direction. Finally it performs a randomn
   *  rotation about the z-axis.
   */
  inline LorentzRotation rotateToZ(Lorentz5Momentum);

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
  bool _useMEforT2;

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
  typedef Herwig::QTildeMECorrection NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the TopDecayMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::TopDecayMECorrection>
  : public ClassTraitsBase<Herwig::TopDecayMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::TopDecayMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * TopDecayMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class TopDecayMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "TopDecayMECorrection.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TopDecayMECorrection.tcc"
#endif

#endif /* HERWIG_TopDecayMECorrection_H */
