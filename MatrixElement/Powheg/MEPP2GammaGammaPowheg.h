// -*- C++ -*-
#ifndef HERWIG_MEPP2GammaGammaPowheg_H
#define HERWIG_MEPP2GammaGammaPowheg_H
//
// This is the declaration of the MEPP2GammaGammaPowheg class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/Shower/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2GammaGammaPowheg class.
 *
 * @see \ref MEPP2GammaGammaPowhegInterfaces "The interfaces"
 * defined for MEPP2GammaGammaPowheg.
 */
class MEPP2GammaGammaPowheg: public HwMEBase {

  enum DipoleType {           IIQCD1=2,           IIQCD2=4, 
		    IFQED1=5, FIQED1=6, IFQED2=7, FIQED2=8 };

  enum RadiationType {Subtraction,Hard,Shower};

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEPP2GammaGammaPowheg();
  //@}
 
  /** @name Member functions for the generation of hard QCD radiation */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return ISR;}

  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr,
				      vector<ShowerInteraction::Type>);
  //@}

public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \f$\alpha_S\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \f$\alpha_{EW}\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

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

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
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
   *  Calculate of the full next-to-leading order weight
   */
  virtual double NLOWeight() const;

  /**
   *  Leading-order matrix element for \f$q\bar q\to \gamma\gamma\f$
   */
  double loGammaGammaME(const cPDVector & particles,
			const vector<Lorentz5Momentum> & momenta,
			bool first=false) const;

  /**
   *  Leading-order matrix element for \f$qg\to \gamma q\f$
   */
  double loGammaqME(const cPDVector & particles,
		const vector<Lorentz5Momentum> & momenta,
		bool first=false) const;

  /**
   *  Leading-order matrix element for \f$g\bar q\to \gamma \bar q\f$
   */
  double loGammaqbarME(const cPDVector & particles,
		       const vector<Lorentz5Momentum> & momenta,
		       bool first=false) const;
  
  /**
   *  Leading-order matrix element for \f$q\bar q\to \gamma g\f$
   */
  double loGammagME(const cPDVector & particles,
		    const vector<Lorentz5Momentum> & momenta,
		    bool first=false) const;

  /**
   *  Real emission matrix element for \f$q\bar q\to \gamma \gamma g\f$
   */
  InvEnergy2 realGammaGammagME(const cPDVector & particles,
			       const vector<Lorentz5Momentum> & momenta,
			       DipoleType dipole, RadiationType rad,
			       bool first=false) const;
  
  /**
   *  Real emission matrix element for \f$qg\to \gamma \gamma q\f$
   */
  InvEnergy2 realGammaGammaqME(const cPDVector & particles,
			       const vector<Lorentz5Momentum> & momenta,
			       DipoleType dipole, RadiationType rad,
			       bool first=false) const;
  
  /**
   *  Real emission matrix element for \f$g\bar q\to \gamma \gamma \bar q\f$
   */
  InvEnergy2 realGammaGammaqbarME(const cPDVector & particles,
				  const vector<Lorentz5Momentum> & momenta,
				  DipoleType dipole, RadiationType rad,
				  bool first=false) const;
  
  /**
   *  The dipole subtractedvirtual contribution
   */
  double subtractedVirtual() const;

  /**
   *  Subtracted real contribution
   */
  double subtractedReal(pair<double,double> x, double z,
			double zJac, double oldqPDF, double newqPDF,
			double newgPDF,bool order) const;

  /**
   *  Calculate of the collinear counterterms
   */
  //@{
  /**
   *  Quark collinear counter term
   */
  double collinearQuark(double x, Energy2 mu2, double jac, double z,
			double oldPDF, double newPDF) const;

  /**
   *  Gluon collinear counter term
   */
  double collinearGluon(Energy2 mu2, double jac, double z,
			double oldPDF, double newPDF) const;
  //@}

  /**
   * The real matrix element divided by \f$2 g_S^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  double realME(const cPDVector & particles,
		const vector<Lorentz5Momentum> & momenta) const;

  /**
   *  Generate hard QCD emission
   */
  HardTreePtr hardQCDEmission(vector<ShowerProgenitorPtr>);

  /**
   *  Generate hard QED emission 
   */
  HardTreePtr hardQEDEmission(vector<ShowerProgenitorPtr>);

  /**
   *  The supression function
   */
  pair<double,double> supressionFunction(Energy pT,Energy scale) const {
    if(supressionScale_==0) scale = lambda_;
    Energy2 scale2 = sqr(scale), pT2 = sqr(pT);
    switch( supressionFunction_ ) {
    case 0:
      return make_pair(1.,0.);
    case 1:
      if(pT < scale ) return make_pair(1.,0.);
      else            return make_pair(0.,1.);
    case 2:
      return make_pair(scale2/(pT2+scale2),pT2/(pT2+scale2));
    default:
      assert(false);
    }
  }


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
  MEPP2GammaGammaPowheg & operator=(const MEPP2GammaGammaPowheg &);

private:

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   *   FFGVertex
   */
  AbstractFFVVertexPtr FFGvertex_;
  //@}

  /**
   *  Kinematic variables for the real radiation
   */
  //@{
  /**
   *  First  variable
   */
  mutable double zTilde_;

  /**
   *  Second variable
   */
  mutable double vTilde_;

  /**
   *  Azimuthal angle
   */
  mutable double phi_;
  //@}

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Pointer to the gluon ParticleData object
   */
  tcPDPtr gluon_;

  /**
   *  Processes
   */
  unsigned int process_;

  /**
   *  Processes
   */
  unsigned int threeBodyProcess_;

  /**
   *  Allowed flavours of the incoming quarks
   */
  int maxflavour_;

  /**
   *  Factor for \f$C_F\f$ dependent pieces
   */
  mutable double CFfact_;

  /**
   *  Factor for \f$T_R\f$ dependent pieces
   */
  mutable double TRfact_;

  /**
   *  Strong coupling
   */
  mutable double alphaS_;

  /**
   *  Use a fixed value of \f$\alpha_S\f$
   */
  bool fixedAlphaS_;

  /**
   *  Electromagnetic coupling
   */
  mutable double alphaEM_;

  /**
   *  Leading-order matrix element
   */
  mutable double loME_;

  /**
   *  The matrix element
   */
  mutable ProductionMatrixElement me_;

  /**
   *  the selected dipole
   */
  mutable DipoleType dipole_;

  /**
   *  Supression Function
   */
  //@{
  /**
   *  Choice of the supression function
   */
  unsigned int supressionFunction_;

  /**
   *  Choice for the scale in the supression function
   */
  unsigned int supressionScale_;

  /**
   *  Scalar for the supression function
   */
  Energy lambda_;
  //@}


  /**
   *  Hard emission stuff
   */
  //@{
  /**
   *  Whether the quark is in the + or - z direction
   */
  bool quarkplus_;
  //@}

  /**
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> beams_;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> partons_;
  //@}

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preQCDqqbarq_;
  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preQCDqqbarqbar_;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double preQCDqg_;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double preQCDgqbar_;

  double preQEDqqbarq_;
  double preQEDqqbarqbar_;
  double preQEDqgq_;
  double preQEDgqbarqbar_;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> prefactor_;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy minpT_;

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaQCD_;

  /**
   *  Pointer to the object calculating the EM
   */
  ShowerAlphaPtr alphaQED_;

  /**
   *  Scale choice
   */
  unsigned int scaleChoice_;

  /**
   *  Scale factor
   */
  double scalePreFactor_;

};

}

#endif /* HERWIG_MEPP2GammaGammaPowheg_H */
