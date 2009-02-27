// -*- C++ -*-
#ifndef HERWIG_DISMECorrection_H
#define HERWIG_DISMECorrection_H
//
// This is the declaration of the DISMECorrection class.
//

#include "QTildeMECorrection.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DISMECorrection class implements the matrix element correction for DIS events.
 *
 * @see \ref DISMECorrectionInterfaces "The interfaces"
 * defined for DISMECorrection.
 */
class DISMECorrection: public QTildeMECorrection {

/**
 *  Typedef for BeamParticleData pointers
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  inline DISMECorrection();

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
			 double & final,EvolverPtr evolver);

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

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
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

private:

  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  inline double generateComptonPoint(double &xp, double & zp);

  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  inline double generateBGFPoint(double &xp, double & zp);

  /**
   *  Calculate the coefficient A for the correlations
   */
  inline double A(tcPDPtr qin, tcPDPtr qout, tcPDPtr lin, tcPDPtr lout);

  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param A \f$\mathcal{A}\f$
   * @param l \f$l=2/y_B-1\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  inline vector<double> ComptonME(double xp, double x2, double xperp,
				  double A, double l, bool norm);

  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_3\f$
   * @param x3 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param A \f$\mathcal{A}\f$
   * @param l \f$l=2/y_B-1\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  inline vector<double> BGFME(double xp, double x2, double x3, double xperp,
			      double A, double l, bool norm);


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<DISMECorrection> initDISMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISMECorrection & operator=(const DISMECorrection &);

private:

  /**
   *   Parameter to control matrix element used
   */
  bool _meopt;

  /**
   *  Radiation enhancement factors
   */
  //@{
  /**
   *  Enchancement factor for ISR
   */
  double _initial;

  /**
   *  Enchancement factor for FSR
   */
  double _final;
  //@}

  /**
   *   Parameters for the phase-space sampling
   */
  //@{
  /**
   *   Relative fraction of compton and BGF processes to generate
   */
  double _procprob;

  /**
   *  Integral for compton process
   */
  double _comptonint;

  /**
   *  Integral for BGF process
   */
  double _bgfint;
  //@}

  /**
   *  Electroweak parameters
   */
  //@{
  /**
   *  \f$\sin\theta_W\f$
   */
  double _sinW;

  /**
   *  \f$\cos\theta_W\f$
   */
  double _cosW;

  /**
   *  The square of the Z mass
   */
  Energy2 _mz2;

  /**
   *  The coefficient for the correlations
   */
  double _acoeff;
  //@}

  /**
   *  Parameters for the point being generated
   */
  //@{
  /**
   *   \f$Q^2\f$
   */
  Energy2 _q2;

  /**
   *  
   */
  double _l;
  //@}

  /**
   *  Testing of weights etc
   */
  //@{
    /**
   *  Number of weights greater than 1
   */
  unsigned int _nover;

  /**
   *  Number of attempts
   */
  unsigned int _ntry;

  /**
   *  Number which suceed
   */
  unsigned int _ngen;

  /**
   *  Maximum weight
   */
  pair<double,double> _maxwgt;

  /**
   *   points for the compton process
   */
  vector<pair<double,double> > _compton,_comptonover;

  /**
   *   points for the BGF process
   */
  vector<pair<double,double> > _bgf,_bgfover;

  /**
   *  Analysis of the x_B dependence
   */
  vector<pair<double,double> > _comptonxb,_bgfxb;
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DISMECorrection. */
template <>
struct BaseClassTrait<Herwig::DISMECorrection,1> {
  /** Typedef of the first base class of DISMECorrection. */
  typedef Herwig::QTildeMECorrection NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DISMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DISMECorrection>
  : public ClassTraitsBase<Herwig::DISMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DISMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DISMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class DISMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "DISMECorrection.icc"

#endif /* HERWIG_DISMECorrection_H */
