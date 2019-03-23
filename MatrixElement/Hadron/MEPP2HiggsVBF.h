// -*- C++ -*-
#ifndef HERWIG_MEPP2HiggsVBF_H
#define HERWIG_MEPP2HiggsVBF_H
//
// This is the declaration of the MEPP2HiggsVBF class.
//

#include "Herwig/MatrixElement/MEfftoffH.h"
#include "Herwig/Shower/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2HiggsVBF class provides the matrix elements for the
 * production of the Higgs boson via the vector boson fusion mechanism
 * in hadron collisions
 *
 * @see \ref MEPP2HiggsVBFInterfaces "The interfaces"
 * defined for MEPP2HiggsVBF.
 */
class MEPP2HiggsVBF: public MEfftoffH {

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
    PPtr incoming;
    
    /**
     *  The outgoing particle
     */
    PPtr outgoing;
    
    /**
     *  The PDF
     */
    tcPDFPtr pdf;
  };

public:

  /**
   * The default constructor.
   */
  MEPP2HiggsVBF();

  /** @name Virtual functions required by the MEBase class. */
  //@{

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;
  //@}

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return Both;}

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

  /**
   *  Apply the POWHEG style correction
   */
  virtual RealEmissionProcessPtr generateHardest(RealEmissionProcessPtr,
						 ShowerInteraction);
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

  /**
   *  Generate the hardest emission in the POWHEG approach
   */
  //@{
  /**
   *  Generate a Compton process
   */
  void generateCompton(unsigned int system);

  /**
   *  Generate a BGF process
   */
  void generateBGF(unsigned int system);

  /**
   *  Matrix element piece for the Compton process
   */
  double comptonME(unsigned int system,
		   double xT,double xp, double zp, double phi);

  /**
   *  Matrix element piece for the Compton process
   */
  double BGFME(unsigned int system, 
	       double xT,double xp, double zp, double phi);
  
  /**
   *  Leading order matrix element
   */
  Energy4 loMatrixElement(const Lorentz5Momentum &p1,
			  const Lorentz5Momentum &p2,
			  const Lorentz5Momentum &q1,
			  const Lorentz5Momentum &q2,
			  double G1, double G2) const;
  //@}

  /**
   *  Generate the hard emission in the old-fashioned matrix element correction approach
   */
  //@{
  /**
   * Generate the values of \f$x_p\f$ and \f$z_p\f$
   * @param xp The value of xp, output
   * @param zp The value of zp, output
   */
  double generateComptonPoint(double &xp, double & zp);

  /**
   * Generate the values of \f$x_p\f$ and \f$z_p\f$
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
   * @param l Scaled momentum of incoming spectator
   * @param m Scaled momentum of outgoing spectator
   *
   */
  vector<double> ComptonME(double xp, double x2, double xperp,
			   LorentzVector<double> l,
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
   * @param l Scaled momentum of incoming spectator
   * @param m Scaled momentum of outgoing spectator
   *
   */
  vector<double> BGFME(double xp, double x2, double x3, double xperp,
		       LorentzVector<double> l,
		       LorentzVector<double> m);
  
  /**
   *  Calculate the coefficient A for the correlations
   */
  double A(tcPDPtr qin1, tcPDPtr qout1, tcPDPtr qin2, tcPDPtr qout2);
  //@}

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
  virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsVBF & operator=(const MEPP2HiggsVBF &) = delete;

private:

  /**
   *  Parameters for the hard POWHEG emission
   */
  //@{
  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alpha_;

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

  /**
   *  Gluon particle data object
   */
  PDPtr gluon_;
  //@}

  /**
   *  Properties of the emission
   */
  //@{
  /**
   *  Beam particle
   */
  tcBeamPtr beam_[2];

  /**
   *  PDF object
   */
  tcPDFPtr pdf_[2];

  /**
   *  Partons
   */
  tcPDPtr partons_[2][4];

  /**
   *  q
   */
  Lorentz5Momentum q_[2];

  /**
   *  \f$Q^2\f$
   */
  Energy2 q2_[2];

  /**
   *  Coupling factor
   */
  double acoeff_;

  /**
   *  Lorentz vectors for the matrix element
   */
  LorentzVector<double> l_;

  /**
   *  Lorentz vectors for the matrix element
   */
  LorentzVector<double> m_;

  /**
   *  Born momentum fraction
   */
  double xB_[2];

  /**
   *  Rotation to the Breit frame
   */
  LorentzRotation rot_[2];

  /**
   *  Quark momenta for spectator system
   */
  Lorentz5Momentum pother_[2][2];

  /**
   *  Quark momenta for emitting system
   */
  Lorentz5Momentum psystem_[2][2];

  /**
   *  Higgs momenta
   */
  Lorentz5Momentum phiggs_[2];

  /**
   *  Transverse momenta for the compton emissions
   */
  Energy pTCompton_[2];

  /**
   *  Transverse momenta for the BGF emissions
   */
  Energy pTBGF_[2];

  /**
   *  Whether the Compton radiation is ISR or FSR
   */
  bool ComptonISFS_[2];

  /**
   *  Momenta of the particles for a compton emission
   */
  vector<Lorentz5Momentum> ComptonMomenta_[2];
  
  /**
   *  Momenta of the particles for a BGF emission
   */
  vector<Lorentz5Momentum> BGFMomenta_[2];

  /**
   *  the systems
   */
  vector<tChannelPair> systems_;

  /**
   *  Higgs boson
   */
  PPtr higgs_;
  //@}

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

  /**
   *  Number of weights greater than 1
   */
  unsigned int nover_;
  
  /**
   *  Maximum weight
   */
  pair<double,double> maxwgt_;
  //@}

};

}

#endif /* HERWIG_MEPP2HiggsVBF_H */
