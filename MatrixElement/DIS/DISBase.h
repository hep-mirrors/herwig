// -*- C++ -*-
#ifndef HERWIG_DISBase_H
#define HERWIG_DISBase_H
//
// This is the declaration of the DISBase class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Shower/ShowerAlpha.h"

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
   * The default constructor.
   */
  virtual ~DISBase();

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
  virtual void initializeMECorrection(RealEmissionProcessPtr, double &,
				      double & );

  /**
   *  Apply the hard matrix element correction to a given hard process or decay
   */
  virtual RealEmissionProcessPtr applyHardMatrixElementCorrection(RealEmissionProcessPtr);

  /**
   * Apply the soft matrix element correction
   * @param parent The initial particle in the current branching
   * @param progenitor The progenitor particle of the jet
   * @param fs Whether the emission is initial or final-state
   * @param highestpT The highest pT so far in the shower
   * @param ids ids of the particles produced in the branching
   * @param z The momentum fraction of the branching
   * @param scale the evolution scale of the branching
   * @param pT The transverse momentum of the branching
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(PPtr parent,
				     PPtr progenitor,
				     const bool & fs,
				     const Energy & highestpT,
				     const vector<tcPDPtr> & ids,
				     const double & z,
				     const Energy & scale,
				     const Energy & pT);
  //@}

  /**
   *  Members for the POWHEG stype correction
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return Both;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr,
						 ShowerInteraction);
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;
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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISBase & operator=(const DISBase &) = delete;

protected:

  /**
   *  The NLO weight
   */
  double NLOWeight() const;

  /**
   *  Calculate the coefficient A for the correlations
   */
  virtual double A(tcPDPtr lin, tcPDPtr lout, tcPDPtr qin, tcPDPtr qout,
		   Energy2 scale) const =0;

  /**
   *  Members for the matrix element correction
   */
  //@{
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
   * @param norm Normalise to the large $l$ value of the ME
   */
  vector<double> ComptonME(double xp, double x2, double xperp,
			   bool norm);
  
  /**
   *  Return the coefficients for the matrix element piece for
   *  the QCD compton case. The output is the \f$a_i\f$ coefficients to 
   *  give the function as 
   *  \f$a_0+a_1\cos\phi+a_2\sin\phi+a_3\cos^2\phi+a_4\sin^2\phi\f$
   * @param xp \f$x_p\f$
   * @param x2 \f$x_3\f$
   * @param x3 \f$x_2\f$
   * @param xperp \f$x_\perp\f$
   * @param norm Normalise to the large $l$ value of the ME
   */
  vector<double> BGFME(double xp, double x2, double x3, double xperp,
		       bool norm);
  //@}

  /**
   *  Members for the POWHEG correction
   */
  //@{
  /**
   *  Generate a Compton process
   */
  void generateCompton();

  /**
   *  Generate a BGF process
   */
  void generateBGF();
  //@}

private:

  /**
   *  Parameters for the matrix element correction
   */
  //@{
  /**
   *  Enchancement factor for ISR
   */
  double initial_;

  /**
   *  Enchancement factor for FSR
   */
  double final_;

  /**
   *   Relative fraction of compton and BGF processes to generate
   */
  double procProb_;

  /**
   *  Integral for compton process
   */
  double comptonInt_;

  /**
   *  Integral for BGF process
   */
  double bgfInt_;
  //@}

  /**
   *  Parameters for the POWHEG correction
   */
  //@{
  /**
   *  Weight for the compton channel
   */
  double comptonWeight_;

  /**
   *  Weight for the BGF channel
   */
  double BGFWeight_;

  /**
   *  Minimum value of \f$p_T\f$
   */
  Energy pTmin_;
  //@}

  /**
   *  Parameters for the point being generated
   */
  //@{
  /**
   *   \f$Q^2\f$
   */
  Energy2 q2_;

  /**
   *  
   */
  double l_;

  /**
   *  Borm momentum fraction
   */
  double xB_;

  /**
   *  Beam particle
   */
  tcBeamPtr beam_;

  /**
   *  Partons
   */
  tcPDPtr partons_[2];

  /**
   *  Leptons
   */
  tcPDPtr leptons_[2];

  /**
   *  PDF object
   */
  tcPDFPtr pdf_;
  /**
   *  Rotation to the Breit frame
   */
  LorentzRotation rot_;

  /**
   *  Lepton momenta
   */
  Lorentz5Momentum pl_[2];

  /**
   *  Quark momenta
   */
  Lorentz5Momentum pq_[2];

  /**
   *  q
   */
  Lorentz5Momentum q_;

  /**
   *  Compton parameters
   */
  Energy pTCompton_;
  bool ComptonISFS_;
  vector<Lorentz5Momentum> ComptonMomenta_;

  /**
   *  BGF parameters
   */
  Energy pTBGF_;
  vector<Lorentz5Momentum> BGFMomenta_;
  //@}

  /**
   *  The coefficient for the correlations
   */
  double acoeff_;

  /**
   *  Coupling
   */
  ShowerAlphaPtr alpha_;

  /**
   *  Gluon particle data object
   */
  PDPtr gluon_;

private:

  /**
   *  The radiative variables
   */
  //@{
  /**
   *  The \f$x_p\f$ or \f$z\f$ real integration variable
   */
  double xp_;
  //@}

  /**
   *  The hadron
   */
  tcBeamPtr hadron_;

  /**
   * Selects a dynamic or fixed factorization scale
   */
  unsigned int scaleOpt_;

  /**
   * The factorization scale 
   */
  Energy muF_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Jacobian for \f$x_p\f$ integral
   */
  double jac_;

};

}

#endif /* HERWIG_DISBase_H */
