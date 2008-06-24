// -*- C++ -*-
#ifndef HERWIG_VBFMECorrection_H
#define HERWIG_VBFMECorrection_H
//
// This is the declaration of the VBFMECorrection class.
//

#include "QTildeMECorrection.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "VBFMECorrection.fh"

namespace Herwig {

using namespace ThePEG;


/**
 *  Typedef for BeamParticleData pointers
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

/**
 *  Struct to contain the hadronic system 
 */
struct tChannelPair{

  /**
   *  The hadron
   */
  PPtr hadron;

  /**
   *  The beam particle data object
   */
  tcBeamPtr beam;

  /**
   *  The incoming particle
   */
  ShowerParticlePtr incoming;

  /**
   *  The outgoing particle
   */
  ShowerParticlePtr outgoing;

  /**
   *  The PDF
   */
  tcPDFPtr pdf;
};

/**
 * The VBFMECorrection class implements the matrix element correction
 * for VBF processes
 *
 * @see \ref VBFMECorrectionInterfaces "The interfaces"
 * defined for VBFMECorrection.
 */
class VBFMECorrection: public QTildeMECorrection {

public:

  /**
   * The default constructor.
   */
  inline VBFMECorrection();

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

  /**
   *  Calculate the coefficient A for the correlations
   */
  inline double A(tcPDPtr qin1, tcPDPtr qout1, tcPDPtr qin2, tcPDPtr qout2);


  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param A \f$\mathcal{A}\f$
   * @param l Scaled momentum of incoming spectator
   * @param m Scaled momentum of outgoing spectator
   *
   */
  inline vector<double> ComptonME(double xp, double x2, double xperp,
				  double A, LorentzVector<double> l,
				  LorentzVector<double> m);

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
   * @param l Scaled momentum of incoming spectator
   * @param m Scaled momentum of outgoing spectator
   *
   */
  inline vector<double> BGFME(double xp, double x2, double x3, double xperp,
			      double A,  LorentzVector<double> l,
			      LorentzVector<double> m);

  /**
   * Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  inline double generateComptonPoint(double &xp, double & zp);

  /**
   * Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  inline double generateBGFPoint(double &xp, double & zp);


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
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<VBFMECorrection> initVBFMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  VBFMECorrection & operator=(const VBFMECorrection &);

private:

  /**
   *   Relative fraction of compton and BGF processes to generate
   */
  double _procprob;

  /**
   *  Integral for compton process
   */
  double _comptonint;

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
  //@}
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of VBFMECorrection. */
template <>
struct BaseClassTrait<Herwig::VBFMECorrection,1> {
  /** Typedef of the first base class of VBFMECorrection. */
  typedef Herwig::QTildeMECorrection NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the VBFMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::VBFMECorrection>
  : public ClassTraitsBase<Herwig::VBFMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::VBFMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * VBFMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class VBFMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMPI.so HwMPIPDF.so HwRemDecayer.so HwShower.so"; }
};

/** @endcond */

}

#include "VBFMECorrection.icc"

#endif /* HERWIG_VBFMECorrection_H */
