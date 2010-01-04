// -*- C++ -*-
#ifndef HERWIG_DISBase_H
#define HERWIG_DISBase_H
//
// This is the declaration of the DISBase class.
//

#include "Herwig++/MatrixElement/HwMEBase.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The DISBase class is the base class for the implementation
 * of DIS type processes including corrections in both the old
 * fashioned matrix element and POWHEG approaches
 *
 * @see \ref DISBaseInterfaces "The interfaces"
 * defined for DISBase.
 */
class DISBase: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  DISBase();

  /**
   *  Members for the old-fashioned matrix element correction
   */
  //@{
  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return true;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(ShowerTreePtr, double &,
				      double & );

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
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr,
				     ShowerParticlePtr,Branching);
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
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<DISBase> initDISBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISBase & operator=(const DISBase &);

protected:

  /**
   *  Calculate the coefficient A for the correlations
   */
  virtual double A(tcPDPtr qin, tcPDPtr qout, tcPDPtr lin, tcPDPtr lout,
		   Energy2 scale) =0;

  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateComptonPoint(double &xp, double & zp);

  /**
   *  Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateBGFPoint(double &xp, double & zp);

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
  vector<double> ComptonME(double xp, double x2, double xperp,
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
  vector<double> BGFME(double xp, double x2, double x3, double xperp,
		       double A, double l, bool norm);

private:

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
   *  The coefficient for the correlations
   */
  double _acoeff;

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
   *  Coupling
   */
  ShowerAlphaPtr _alpha;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DISBase. */
template <>
struct BaseClassTrait<Herwig::DISBase,1> {
  /** Typedef of the first base class of DISBase. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DISBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::DISBase>
  : public ClassTraitsBase<Herwig::DISBase> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::DISBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MENeutralCurrentDIS is implemented. It may also include several, space-separated,
   * libraries if the class MENeutralCurrentDIS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEDIS.so"; }
};

/** @endcond */

}

#endif /* HERWIG_DISBase_H */
