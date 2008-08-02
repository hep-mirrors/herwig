// -*- C++ -*-
//
// MEPP2HiggsPowheg.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2HiggsPowheg_H
#define HERWIG_MEPP2HiggsPowheg_H
//
// This is the declaration of the MEPP2HiggsPowheg class.
//

#include "ThePEG/MatrixElement/MEBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "Herwig++/PDT/SMHiggsMassGenerator.h"
#include "Herwig++/Models/StandardModel/StandardModel.h"
#include "Herwig++/MatrixElement/General/ProductionMatrixElement.h"
#include "ThePEG/PDF/BeamParticleData.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2HiggsPowheg class implements the matrix element for the process
 * pp->Higgs with different Higgs shape prescriptions (see details in hep-ph/9505211)
 * and the NLL corrected Higgs width (see details in the FORTRAN HERWIG manual).
 *
 * @see \ref MEPP2HiggsPowhegInterfaces "The interfaces"
 * defined for MEPP2HiggsPowheg.
 */
class MEPP2HiggsPowheg: public MEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2HiggsPowheg();

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
  virtual CrossSection dSigHatDR() const;

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
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \f$ ]0,1[ \f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

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
  virtual Selector<const ColourLines *> colourGeometries(tcDiagPtr diag) const;

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
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
   *  Parameters for the NLO weight
   */

  /**
   *  The colour factors
   */
  double CF_ , CA_ , TR_;

  /**
   * Number of light flavours (in the beta function beta0_)
   */
  double nlf_;

  /**
   * (Proportional to) The beta function
   */
  double beta0_;

  //@{
  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Whether to use a fixed or a running QCD coupling for the NLO weight
   */
  unsigned int nlo_alphaS_opt_;

  /**
   *  The value of alphaS to use for the nlo weight if nlo_alphaS_opt_=1
   */
  double fixed_alphaS_;
  //@}

  /**
   *  Radiation variables
   */
  //@{
  /**
   *   The \f$\tilde{x}\f$ variable
   */
  double xt_;

  /**
   *  The \f$y\f$ angular variable
   */
  double y_;
  //@}

  /**
   *  Values of the PDF's before radiation
   */
  mutable double lo_lumi_;

  /**
   *  The value of the leading order gg->H matrix element
   */
  mutable double lo_ggME_;

  /**
   * Function to set the born variables. 
   */
  void get_born_variables() const    ;

  /**
   * The invariant mass of the lo final state. 
   */
  mutable Energy2 p2_     , s2_      ;

  /**
   * The squared masses of the lo final state particles p1 and p2. 
   */
  mutable Energy2 p12_    , p22_     ;

  /**
   * The polar and azimuthal angles respectively defining a two body lo event. 
   */
  mutable double  theta1_ , theta2_  ;

  /**
   *  The momentum fraction of the plus and minus partons in the Born process
   */
  mutable double xbp_, xbm_;

  /**
   *  The sqrt(1-xbp_) and sqrt(1-xbm_) respectively
   */
  mutable double etabarp_, etabarm_;

  /**
   *  The ParticleData object for the plus and minus lo partons
   */
  mutable tcPDPtr a_lo_, b_lo_;

  /**
   *  The BeamParticleData object for the plus and minus direction hadrons
   */
  mutable Ptr<BeamParticleData>::transient_const_pointer hadron_A_;
  mutable Ptr<BeamParticleData>::transient_const_pointer hadron_B_;

  /**
   *  The value of \f$\alpha_S\f$ used for the calculation
   */
  mutable double alphaS_;

  /**
   * Calculate the correction weight with which leading-order
   * configurations are re-weighted.
   */
  double NLOweight() const;

  /**
   * Invariants required for the evaluation of next-to-leading order
   * quantities (Frixione et al. NPB.383 WZ production at colliders). 
   */
  Energy2 s(double xt, double y)      const ; 
  Energy2 tk(double xt, double y)     const ;
  Energy2 uk(double xt, double y)     const ; 
  double  betax(double xt, double y)  const ; 
  double  v1(double xt, double y)     const ; 
  double  v2(double xt, double y)     const ; 
  double  cpsi(double xt, double y)   const ; 
  double  cpsipr(double xt, double y) const ; 
  Energy2 q1(double xt, double y)     const ;
  Energy2 q2(double xt, double y)     const ;
  Energy2 q1hat(double xt, double y)  const ; 
  Energy2 q2hat(double xt, double y)  const ; 
  Energy2 w1(double xt, double y)     const ; 
  Energy2 w2(double xt, double y)     const ;  

  /**
   * Calculate the minimum of \f$x\f$. 
   */
  double xbar(double y) const;

  /**
   * Calculate auxiliary function of \f$\bar{x}(y)\f$, \f$\bar{\eta}(y)\f$. 
   */
  double etabar(double y) const;

  /**
   * Calculate the variable \f$x=p^{2}/s\f$ from the integration variables. 
   */
  double x(double xt, double y) const;

  /**
   * Calculate the momentum fraction of the plus and minus partons. 
   */
  double xp(double x, double y) const;
  double xm(double x, double y) const;

  /**
   * Calculate the ratio of the NLO luminosity to the LO
   * luminosity function for the \f$q\bar{q}\f$ initiated channel. 
   */
  double Lhat_ab(tcPDPtr a, tcPDPtr b, double x, double y) const;

  /**
   * Calculate the universal soft-virtual contribution to the NLO weight. 
   */
  double Vtilde_universal() const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_qq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$gg\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_qg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Ctilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

  /**
   * The regular part of the virtual correction matrix element(s) 
   */
  double M_V_regular() const;

  /**
   * The matrix element q + qbar -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qqbar(double xt, double y) const;

  /**
   * The matrix element qbar + q -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_qbarq(double xt, double y) const;

  /**
   * The matrix element g + g    -> n + g times tk*uk 
   */
  Energy2 t_u_M_R_gg(double xt, double y) const;

  /**
   * The matrix element q + g    -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_qg(double xt, double y) const;

  /**
   * The matrix element g + q    -> n + q times tk*uk 
   */
  Energy2 t_u_M_R_gq(double xt, double y) const;

  /**
   * Function for calculation of the \f$q\bar{q}\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qqbar_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;
  /**
   * Function for calculation of the \f$\bar{q}q\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qbarq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;
  /**
   * Function for calculation of the \f$qq\f$ 
   * initiated real contribution.
   */
  double Rtilde_Ltilde_gg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;
  /**
   * Function for calculation of the \f$qg\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_qg_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;
  /**
   * Function for calculation of the \f$gq\f$ initiated real
   * contribution.
   */
  double Rtilde_Ltilde_gq_on_x(tcPDPtr a,tcPDPtr b,double xt,double y) const;

protected:
  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const { return new_ptr(*this); }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const { return new_ptr(*this); }
  //@}

protected:
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit() throw(InitException);

private:
  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2HiggsPowheg> initMEPP2HiggsPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsPowheg & operator=(const MEPP2HiggsPowheg &);
  //@}

private:
  /**
   *  Members to return the matrix elements for the different subprocesses
   */
  //@{

  /**
   * Calculates the matrix element for the process g,g->h (via quark loops)
   * @param g1 a vector of wave functions of the first incoming gluon
   * @param g2 a vector of wave functions of the second incoming gluon
   * @param calc Whether or not to calculate the matrix element for spin correlations
   * @return the amlitude value.
   */
  double ggME(vector<VectorWaveFunction> g1,
              vector<VectorWaveFunction> g2,
              ScalarWaveFunction &, 
              bool calc) const;

  /**
   * Calculates the matrix element for the process q,qbar->h
   * @param fin a vector of quark spinors
   * @param ain a vector of anti-quark spinors
   * @param calc Whether or not to calculate the matrix element for spin correlations
   * @return the amlitude value.
   */
  double qqME(vector<SpinorWaveFunction> & fin, 
              vector<SpinorBarWaveFunction> & ain, 
              ScalarWaveFunction &, 
              bool calc) const;
  //@}

private:

  /**
   * Selects a dynamic (sHat) or fixed factorization scale
   */
  unsigned int scaleopt_;

  /**
   * The factorization and renormalization scale respectively
   */
  Energy mu_F_, mu_UV_;

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   * Defines the Higgs resonance shape
   */
  unsigned int shapeopt_;

  /**
   * The processes to be included (GG->H and/or qq->H)
   */
  unsigned int processopt_;

  /**
   * Minimum flavour of incoming quarks
   */
  unsigned int minflavouropt_;

  /**
   * Maximum flavour of incoming quarks
   */
  unsigned int maxflavouropt_;

  /**
   * Storage of the diagram weights for the \f$gg\to Hg\f$ subprocess
   */
  mutable double diagwgt[3];

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement me_;

  /**
   * Pointer to the H->2gluons vertex (used in gg->H)
   */
  AbstractVVSVertexPtr hggvertex;

  /**
   * Pointer to the fermion-fermion Higgs vertex (used in qq->H)
   */
  AbstractFFSVertexPtr ffhvertex;

  /**
   * Pointer to the Standard Model instance used in the class
   */
  tcHwSMPtr theSM;

  /**
   *  The mass generator for the Higgs
   */
  SMHiggsMassGeneratorPtr hmass_;

  /**
   *  On-shell mass for the higgs
   */
  Energy mh_;

  /**
   *  On-shell width for the higgs
   */
  Energy wh_;

  /**
   *  Total branching for the allowed decays; this compensates
   *  for the fact that the SMHiggsWidthGenerator currently does
   *  take account of whether decay modes are switched On / Off 
   *  when using the WidthScheme=Fixed option.
   */
  double h_br_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes of MEPP2HiggsPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2HiggsPowheg,1> {
  /** Typedef of the first base class of MEPP2HiggsPowheg. */
  typedef MEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2HiggsPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2HiggsPowheg>
  : public ClassTraitsBase<Herwig::MEPP2HiggsPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2HiggsPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2HiggsPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2HiggsPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegME.so"; }
};

/** @endcond */

}

#include "MEPP2HiggsPowheg.icc"

#endif /* HERWIG_MEPP2HiggsPowheg_H */
