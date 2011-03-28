// -*- C++ -*-
#ifndef HERWIG_MEPP2VGammaPowheg_H
#define HERWIG_MEPP2VGammaPowheg_H
//
// This is the declaration of the MEPP2VGammaPowheg class.
//

#include "Herwig++/MatrixElement/Hadron/MEPP2VGamma.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "Herwig++/Utilities/GSLIntegrator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"  //from VGammaHardGenerator
#include "Herwig++/Shower/Base/ShowerProgenitor.h"  //from VGammaHardGenerator


namespace Herwig {
  
using namespace ThePEG;

class MEPP2VGammaPowheg;
ThePEG_DECLARE_POINTERS(MEPP2VGammaPowheg,VGamPowhegPtr);

/**
 * The MEPP2VGammaPowheg class implements the next-to-leading
 * order matrix elements for $q\bar q \to W^\pm/Z^0\gamma\f$
 * in the Powheg scheme.
 *
 * @see \ref MEPP2VGammaPowhegInterfaces "The interfaces"
 * defined for MEPP2VGammaPowheg.
 */
class MEPP2VGammaPowheg: public MEPP2VGamma {

   /**
   * Typedef for the BeamParticleData object from VGammaHardGenerator
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr; 

  /*
public:
  friend class KP1Integrand;
  friend class KP2Integrand;
  friend class KP3Integrand;  

  */
public:

  /**
   * The default constructor.
   */
  MEPP2VGammaPowheg();

  /** @name Member functions for the generation of hard QCD radiation */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual bool hasPOWHEGCorrection() {return true;}

  /**
   *  Member to generate the Powheg hardest emission from VGammaHardGenerator
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr,vector<ShowerInteraction::Type>);


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

  //parton distribution function of beam partID
  double PDFab(double z, int partID, Energy2 scale, tcBeamPtr hadron_A, tcBeamPtr hadron_B)  const;
  // parton distribution function ratio
  double PDFratio(double z, double x, Energy2 scale, int partID, int fi, tcBeamPtr hadron_A, tcBeamPtr hadron_B)  const;


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

  

  // functionals of PDFs and energy fractions in the collinear remnants:
  double KP1(double zz, Energy2 shat, Energy2 muF2, int partID) const;
  double KP2(double zz, Energy2 shat, Energy2 muF2) const;
  double KP3(double zz, int partID) const;
  double KPpr(double xx, Energy2 shat, Energy2 muF2) const;

  // definitions of H, F^V:
  double Hfunc(Energy2 t, Energy2 s, Energy2 m2) const;
  double FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;
  double FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;
  double FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;

  /**
   *  Generate a \f$q \bar q \to V \gamma \f$ configuration
   */
  void generateQQbarG();

  /**
   *   The ratio of leading to NLO matrix elements for qqbar to gluon V gamma
   */ 
  double QQbarGratio(Lorentz5Momentum k_a, Lorentz5Momentum  k_b, Lorentz5Momentum k_ga, Lorentz5Momentum k_V, 
                     Lorentz5Momentum k_glu, Energy2 scale);

  /**
   * generate qg to q Vgamma events (qgflag=1) and gq_bar to q_bar Vgamma events (qgflag=2)
   */
  void generateQGQ(int qgflag);

  /**
   *   The ratio of leading to NLO matrix elements for gq_bar to q_bar Vgamma
   */
  double QGQratio(Lorentz5Momentum k_a, Lorentz5Momentum  k_b, Lorentz5Momentum k_ga, Lorentz5Momentum k_V, 
		  Lorentz5Momentum k_q, Energy2 scale, int flavorflag);


// Transformation from Born+rad phase space configuration to m+1 phase space of final states
  Lorentz5Momentum InvLortr(Lorentz5Momentum kbar_a, Lorentz5Momentum kbar_b, double xi,
			    Lorentz5Momentum k_i, Lorentz5Momentum kbar_j, int fi) const;

//Transformation from m+1 phase space to Born phase space configuration of final states
  Lorentz5Momentum Lortr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, double xi,
			 Lorentz5Momentum k_i, Lorentz5Momentum k_j, int fi) const;

// momentum of radiated parton transformed from Born phase space
  Lorentz5Momentum radk(Lorentz5Momentum kbar_a, Lorentz5Momentum kbar_b, double xi,
			double vi, double phi, int fi) const;

// momentum of radiated quark transformed from gq->qV Born phase space
  Lorentz5Momentum radkqr(Lorentz5Momentum kbar_ga, Lorentz5Momentum kbar_V, Energy2 MassV2, 
			  double zi, double ui, double phi, double zcut) const;
// Born matrix element
  double MatrBorn(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		  Lorentz5Momentum k_V, double char0, double char1, Energy2 MassV2) const;

// Tree level matrix element for gq->qV
  double MatrQV(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
		Lorentz5Momentum k_V, double charqr, Energy2 MassV2, double alphas, int fi) const;
// Gluon Real radiation matrix element
  double MatrRealGluon(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		       Lorentz5Momentum k_V, Lorentz5Momentum k_glu, Energy2 scale) const;

// Quark Real radiation matrix element
  double MatrRealQuark(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		       Lorentz5Momentum k_V, Lorentz5Momentum k_q, Energy2 scale, int flavorflag) const;

// dipole of gluon radiation
  double Dipole(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		Lorentz5Momentum k_V,double char0,double char1,Energy2 MassV2,Energy2 kig, double alphas,double xi) const;

// dipole of quark radiation in the singular region collinear with final state photon
  double Dipolepqr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		   Lorentz5Momentum k_V, double chargr, Energy2 MassV2, Energy2 kig, double zi, double alphas, int fi) const;



// dipole of quark radiation in the singular region collinear with initial state gluon
  double Dipolegluqr(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		     Lorentz5Momentum k_V, double char0, double char1, Energy2 MassV2, Energy2 kig, double alphas, double xi) const;

// Born matrix element for V quark procduction
  double qgME(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
				 Lorentz5Momentum k_V, Energy2 MassV2, bool calc)  const;

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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2VGammaPowheg> initMEPP2VGammaPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VGammaPowheg & operator=(const MEPP2VGammaPowheg &);

protected:

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;


private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaS_;

  mutable double  alphaQCD_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pTmin_;


  /**
   * The \f$T_R\f$ colour factor
   */
  mutable double _TR;

  /**
   *  The \f$C_F\f$ colour factor
   */
  mutable double _CF;

  mutable double Pi;
  
  mutable double gamnum; 
  //Factorization scale
  mutable Energy2 _muF2;

    //kinematics: s~hat:_ss, t~hat:_tt, u~hat:_uu
  mutable Energy2 _ss;
  mutable Energy2 _tt;
  mutable Energy2 _uu;

  //types of initial states:
  mutable tcPDPtr _partona;
  mutable tcPDPtr _partonb;

    //types of final states:
    mutable tcPDPtr _gluon;
    mutable tcPDPtr _photon;
    mutable tcPDPtr _quark; 
    mutable tcPDPtr _quarkpi;
    mutable tcPDPtr _boson;
    mutable int _idboson;
    mutable int _iflagcq;
    
    //momenta of final states:
    mutable Lorentz5Momentum _p_photon;
    mutable Lorentz5Momentum _p_parton;
    mutable Lorentz5Momentum _p_boson;
    
    //momenta of initial states:
    mutable Lorentz5Momentum _p_partona;
    mutable Lorentz5Momentum _p_partonb;
   
  //CKM matrix square
  mutable double ckm;
  // the third componet of isospin of the quark when vector boson is Z0
  mutable double I3quark;
  /**
   *  The BeamParticleData object for the plus direction hadron
   */
  mutable tcBeamPtr _hadron_A;

  /**
   *  The BeamParticleData object for the  minus direction hadron
   */
  mutable tcBeamPtr _hadron_B;

    //momenta fractions _xa, _xb
    mutable double _xa;
    mutable double _xb;

  //
    //alpha_S:
    mutable double _alphas;
  //sin(thete_W)^2
    mutable double sin2w;
  //electroweak alpha
    mutable double alphae;

    //integrated variable _xtil
    mutable double _xtil;

  //Real radiation integrated variables
    mutable double _y1;
    mutable double _y2;
    mutable double _y3;
 
  // For Real radiation part
    mutable double _phi;
    
  /**
   *  Parameters for the NLO weight
   */
  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int _contrib;
  //@}

  /**
   *  Choice of the scale
   */
  //@{
  /**
   *  Type of scale
   */
  unsigned int _scaleopt;

  /**
   *  Fixed scale if used
   */
  Energy _fixedScale;

  /**
   *  Prefactor if variable scale used
   */
  double _scaleFact;

  //photon fraction cut:
  mutable double fraccut;
  //@}

  /**
   * integrator
   */
  GSLIntegrator _integrator;

  /**
   *  Vertices
   */
  //@{
  /**
   *   FFPVertex
   */
  AbstractFFVVertexPtr FFPvertex_;

  /**
   *   FFWVertex
   */
  AbstractFFVVertexPtr FFWvertex_;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr FFZvertex_;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr WWWvertex_;

  /**
   *  FFG Vertex
   */ 
  AbstractFFVVertexPtr FFGvertex_;
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
   *  Whether the quark is in the + or - z direction
   */
  bool quarkplus_;

  /**
   *  Prefactor for the overestimate for \f$q\bar q\to V \gamma\f$
   */
  double qqgFactor_;

  /**
   *  Prefactor for the overestimate for q g to V gamma q
   */
  double qgFactor_g;
  double qgFactor_p;

  /**
   *  Prefactor for the overestimate for qbar g to V gamma qbar
   */
  double gqbarFactor_g;
  double gqbarFactor_p;

  /**
   * The power, \f$n\f$, for the sampling of initial state partons singular regions
   */
  double power_;

  /**
   * The power, \f$n\f$, for the sampling of final state photon singular region
   */
  double power_photon;

  /**
   *  Born variables
   */
  //@{
  /**
   *  Rapidity of the photon
   */
  double photonRapidity_;

  /**
   *  Rapidity of the gauge boson
   */
  double bosonRapidity_;

  /**
   *  \f$p_T\f$ of the photon
   */
  Energy photonpT_;

  /**
   *  Azimuth of the photon
   */
  double photonAzimuth_;

  /**
   * gauge boson mass square
   */
  mutable Energy2 MV2;

  /**
   * Mass of the boson/photon system
   */
  Energy systemMass_;

  /**
   *  Momentum fractions for the LO process
   */ 
  double x_[2];


  // The charges of the initial state quarks
  mutable double charge0;
  mutable double charge1;
  // The charges of the radiated quarks
  mutable double chargeqr;
  //@}

  /**
   *  Momenta etc for the \f$q\bar q \to V \gamma g \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVqqbar_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammaqqbar_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGqqbar_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQqqbar_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQbarqqbar_;

  /**
   *  The transverse momentum
   */
  Energy pTqqbar_;
  //@}

  /**
   *  Momenta etc for the \f$qg  \to V \gamma q \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVqg_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammaqg_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGqg_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQinqg_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQoutqg_;

  /**
   *  The transverse momentum
   */
  Energy pTqg_;
  //@}

  /**
   *  Momenta etc for the \f$g\bar q \to V \gamma \bar q \f$ process
   */
  //@{
  /**
   *  Momentum of the vector boson
   */
  Lorentz5Momentum pVgqbar_;

  /**
   *  Momentum of the photon
   */
  Lorentz5Momentum pGammagqbar_;

  /**
   *  Momentum of the gluon
   */
  Lorentz5Momentum pGgqbar_;

  /**
   *  Momentum of the incoming quark
   */
  Lorentz5Momentum pQingqbar_;

  /**
   *  Momentum of the incoming antiquark
   */
  Lorentz5Momentum pQoutgqbar_;

  /**
   *  The transverse momentum
   */
  Energy pTgqbar_;

  //** the variables that are passed to the function QQbarGratio :
  Energy2 ssHat_;
  Energy  maTV_;
  Energy pTparton_;
  double phiparton_;
  double yparton_;  

  // ISR flag for QGQ generation
  int QGQISR;
  int QGQISR_gqbar;
  int QGQISR_qg;


};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VGammaPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2VGammaPowheg,1> {
  /** Typedef of the first base class of MEPP2VGammaPowheg. */
  typedef Herwig::MEPP2VGamma NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VGammaPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VGammaPowheg>
  : public ClassTraitsBase<Herwig::MEPP2VGammaPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VGammaPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VGammaPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VGammaPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

/*
namespace Herwig {
  // friend class MEPP2VGammaPowheg; 
  //  using Herwig::MEPP2VGammaPowheg::PDFab;
 struct KP1Integrand  {
   KP1Integrand(tcVGamPowhegPtr MEclass, double zz, Energy2 shat, Energy2 muF2, int partID, tcBeamPtr hadron_A, tcBeamPtr hadron_B) 
     : _meclass(MEclass), _zz(zz), _shat(shat), _muF2(muF2), _pID(partID), _hadron_A(hadron_A), _hadron_B(hadron_B) {}

   double operator ()(double x) const {
     double pdfzx, pdfz;

     pdfzx= _meclass->PDFab(_zz/x, _pID, _muF2, _hadron_A, _hadron_B);
     pdfz=  _meclass->PDFab(_zz, _pID, _muF2, _hadron_A, _hadron_B);
     return log(sqr(1.0-x)*_shat/(x*_muF2))*
       (2.0/(1.0-x) -(1.0+sqr(x))/(1.0-x)*pdfzx/pdfz/x) +2.0/(1.0-x)*log(x);
   }
   typedef double ArgType;
   typedef double ValType;
 
   tcVGamPowhegPtr _meclass;
   double _zz;
   Energy2 _shat;
   Energy2 _muF2;
   int _pID;  
   
 };

 struct KP2Integrand  {
   KP2Integrand(Energy2 shat, Energy2 muF2) 
     : _shat(shat), _muF2(muF2) {}

   double operator ()(double x) const {
     return log(sqr(1.0-x)*_shat/_muF2)*(2.0/(1.0-x));
   }
   typedef double ArgType;
   typedef double ValType;
 
   Energy2 _shat;
   Energy2 _muF2;
 };

 struct KP3Integrand  {
   KP3Integrand(tcVGamPowhegPtr MEclass, double zz, int partID, Energy2 muF2, tcBeamPtr hadron_A, tcBeamPtr hadron_B) 
     : _meclass(MEclass), _zz(zz), _pID(partID), _muF2(muF2), _hadron_A(hadron_A), _hadron_B(hadron_B)  {}

   double operator ()(double x) const {
     double pdfzx, pdfz;
     MEPP2VGammaPowheg VGamPow;
     pdfzx= _meclass->PDFab(_zz/x, _pID, _muF2, _hadron_A, _hadron_B);
     pdfz=  _meclass->PDFab(_zz, _pID, _muF2, _hadron_A, _hadron_B);

     return (1.0-x)*pdfzx/pdfz/x;
   }
   typedef double ArgType;
   typedef double ValType;
 
   tcVGamPowhegPtr _meclass;
   double _zz;
   int _pID;  
 };
}
*/


#endif /* HERWIG_MEPP2VGammaPowheg_H */
