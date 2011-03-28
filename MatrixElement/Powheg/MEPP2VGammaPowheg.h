// -*- C++ -*-
#ifndef HERWIG_MEPP2VGammaPowheg_H
#define HERWIG_MEPP2VGammaPowheg_H
//
// This is the declaration of the MEPP2VGammaPowheg class.
//

#include "Herwig++/MatrixElement/Hadron/MEPP2VGamma.h"
#include "ThePEG/PDF/BeamParticleData.h"


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

public:

  /**
   * The default constructor.
   */
  MEPP2VGammaPowheg();

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
  double PDFab(double z, int partID)  const;
  // parton distribution function ratio
  double PDFratio(double z, double x, int partID)  const;

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

  // definitions of H, F^V:
  double Hfunc(Energy2 t, Energy2 s, Energy2 m2) const;
  double FWfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;
  double FZfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;
  double FVfunc(Energy2 t, Energy2 u, Energy2 s, Energy2 m2) const;

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
  Lorentz5Momentum radkqr(Lorentz5Momentum kbar_ga, Lorentz5Momentum kbar_V, Energy2 MV2, 
			  double zi, double ui, double phi) const;
// Born matrix element
  double MatrBorn(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		  Lorentz5Momentum k_V, double char0, double char1, Energy2 MV2) const;

// Tree level matrix element for gq->qV
  double MatrQV(Lorentz5Momentum p_a, Lorentz5Momentum p_b, Lorentz5Momentum k_q,
		Lorentz5Momentum k_V, Energy2 MV2) const;
  // Gluon Real radiation matrix element
  InvEnergy2 MatrRealGluon(Lorentz5Momentum k_a, Lorentz5Momentum k_b,
			   Lorentz5Momentum k_ga,
			   Lorentz5Momentum k_V, Lorentz5Momentum k_glu) const;
  
// Quark Real radiation matrix element
  double MatrRealQuark(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		       Lorentz5Momentum k_V, Lorentz5Momentum k_q) const;

// dipole of gluon radiation
  InvEnergy2 Dipole(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		    Lorentz5Momentum k_V, double char0,double char1,Energy2 MV2,Energy2 kig,double xi) const;

  double test2to3(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		  Lorentz5Momentum k_1, Lorentz5Momentum k_2) const;
  double test2to3am(Lorentz5Momentum k_a, Lorentz5Momentum k_b, Lorentz5Momentum k_ga,
		  Lorentz5Momentum k_1, Lorentz5Momentum k_2, Energy2 MV2) const;

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
   * The \f$T_R\f$ colour factor
   */
  mutable double TR_;

  /**
   *  The \f$C_F\f$ colour factor
   */
  mutable double CF_;
    
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
    mutable int _iflagcq;
    
    //momenta of final states:
    mutable Lorentz5Momentum _p_photon;
    mutable Lorentz5Momentum _p_parton;
    mutable Lorentz5Momentum _p_boson;
    
    //momenta of initial states:
    mutable Lorentz5Momentum _p_partona;
    mutable Lorentz5Momentum _p_partonb;
   
  //CKM matrix square
  mutable double ckm_;

  /**
   *  The BeamParticleData object for the plus direction hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer _hadron_A;

  /**
   *  The BeamParticleData object for the  minus direction hadron
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer _hadron_B;

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
  //@}

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

#endif /* HERWIG_MEPP2VGammaPowheg_H */
