// -*- C++ -*-
//
// MEPP2Higgs.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2Higgs_H
#define HERWIG_MEPP2Higgs_H
//
// This is the declaration of the MEPP2Higgs class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "Herwig/Shower/Core/Couplings/ShowerAlpha.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2Higgs class implements the matrix element for the process
 * pp->Higgs with different Higgs shape prescriptions (see details in hep-ph/9505211)
 * and the NLL corrected Higgs width (see details in the FORTRAN HERWIG manual).
 *
 * @see \ref MEPP2HiggsInterfaces "The interfaces"
 * defined for MEPP2Higgs.
 */
class MEPP2Higgs: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2Higgs();

  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object.
   */
  virtual void setKinematics() {
    HwMEBase::setKinematics();
    mh2_ = sHat();
  }

public:
 
  /** @name Member functions for the generation of hard QCD radiation */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return ISR;}

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
   * @param initial The particle from the hard process which started the 
   * shower
   * @param parent The initial particle in the current branching
   * @param br The branching struct
   * @return If true the emission should be vetoed
   */
  virtual bool softMatrixElementVeto(ShowerProgenitorPtr initial,
				     ShowerParticlePtr parent,
				     Branching br);

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
  
  /**
   *   Members to calculate the real emission matrix elements
   */
  //@{
  /**
   *  The leading-order matrix element for \f$gg\to H\f$
   */
  Energy4 loME() const;

  /**
   *  The matrix element for \f$gg\to H g\f$
   */
  Energy2 ggME(Energy2 s, Energy2 t, Energy2 u);

  /**
   *  The matrix element for \f$qg\to H q\f$
   */
  Energy2 qgME(Energy2 s, Energy2 t, Energy2 u);

  /**
   *  The matrix element for \f$qbarg\to H qbar\f$
   */
  Energy2 qbargME(Energy2 s, Energy2 t, Energy2 u);
  //@}

  /**
   *  Members to calculate the functions for the loop diagrams
   */
  //@{
  /**
   *  The \f$B(s)\f$ function of NBP339 (1990) 38-66
   * @param s The scale
   * @param mf2 The fermion mass squared.
   */
  Complex B(Energy2 s,Energy2 mf2) const;

  /**
   *  The \f$C(s)\f$ function of NBP339 (1990) 38-66
   * @param s The scale
   * @param mf2 The fermion mass squared.
   */
  complex<InvEnergy2> C(Energy2 s,Energy2 mf2) const;

  /**
   *  The \f$C(s)\f$ function of NBP339 (1990) 38-66
   * @param s The \f$s\f$ invariant
   * @param t The \f$t\f$ invariant
   * @param u The \f$u\f$ invariant
   * @param mf2 The fermion mass squared
   */
  complex<InvEnergy4> D(Energy2 s,Energy2 t, Energy2 u,Energy2 mf2) const;

  /**
   * The integral \f$\int\frac{dy}{y-y_0}\log(a-i\epsilon-b y(1-y))\f$
   * from NBP339 (1990) 38-66.
   * @param a  The parameter \f$a\f$.
   * @param b  The parameter \f$b\f$.
   * @param y0 The parameter \f$y_0\f$.
   */
  Complex dIntegral(Energy2 a, Energy2 b, double y0) const;

  /**
   *  The \f$M_{+++}\f$ matrix element of NBP339 (1990) 38-66.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   * @param i Which of the stored values to use for \f$D(u,t)\f$.
   * @param j Which of the stored values to use for \f$D(u,s)\f$.
   * @param k Which of the stored values to use for \f$D(s,t)\f$.
   * @param i1 Which of the stored values to use for \f$C_1(s)\f$.
   * @param j1 Which of the stored values to use for \f$C_1(t)\f$.
   * @param k1 Which of the stored values to use for \f$C_1(u)\f$.
   */
  complex<Energy> me1(Energy2 s,Energy2 t,Energy2 u, Energy2 mf2,
			     unsigned int i,unsigned int j, unsigned int k,
			     unsigned int i1,unsigned int j1, unsigned int k1) const;

  /**
   *  The \f$M_{++-}\f$ matrix element of NBP339 (1990) 38-66.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  complex<Energy> me2(Energy2 s,Energy2 t,Energy2 u, Energy2 mf2) const;

  /**
   *  The \f$F(x)\f$ function for the leading-order result
   */
  Complex F(double x) const;
  //@}

  /**
   *  Method to extract the PDF weight for quark/antiquark
   * initiated processes and select the quark flavour
   */  
  tPDPtr quarkFlavour(tcPDFPtr pdf, Energy2 scale, double x, tcBeamPtr beam, 
		      double & pdfweight, bool anti);

  /**
   * Return the momenta and type of hard matrix element correction
   * @param gluons The original incoming particles.
   * @param beams The BeamParticleData objects
   * @param higgs The original outgoing higgs
   * @param iemit Whether the first (0) or second (1) particle emitted
   * the radiation
   * @param itype The type of radiated particle (0 is gluon, 1 is quark 
   *              and 2 is antiquark)
   * @param pnew The momenta of the new particles
   * @param xnew The new values of the momentuym fractions
   * @param out The ParticleData object for the outgoing parton
   * @return Whether or not the matrix element correction needs to be applied
   */
  bool applyHard(ParticleVector gluons,
		 vector<tcBeamPtr> beams,
		 PPtr higgs,unsigned int & iemit,
		 unsigned int & itype,vector<Lorentz5Momentum> & pnew,
		 pair<double,double> & xnew,
		 tPDPtr & out);

  /**
   *  generates the hardest emission (yj,p)
   * @param pnew The momenta of the new particles
   * @param emissiontype The type of emission, as for getResult
   * @return Whether not an emission was generated
   */
  bool getEvent(vector<Lorentz5Momentum> & pnew,int & emissiontype);

  /**
   * Returns the matrix element for a given type of process,
   * rapidity of the jet \f$y_j\f$ and transverse momentum \f$p_T\f$
   * @param emis_type the type of emission,
   * (0 is \f$gg\to h^0g\f$, 1 is \f$qg\to h^0q\f$ and 2 is \f$g\bar{q}\to h^0\bar{q}\f$)
   * @param pt The transverse momentum of the jet
   * @param yj The rapidity of the jet
   * @param outParton the outgoing parton
   */
  double getResult(int emis_type, Energy pt, double yj,tcPDPtr & outParton);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2Higgs> initMEPP2Higgs;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2Higgs & operator=(const MEPP2Higgs &);
  //@}

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
   * The value associated to the fixed factorization scale option
   */
  Energy mu_F_;

  /**
   * Defines the Higgs resonance shape
   */
  unsigned int shapeOption_;

  /**
   * The processes to be included (GG->H and/or qq->H)
   */
  unsigned int processOption_;

  /**
   * Minimum flavour of incoming quarks
   */
  int minFlavour_;

  /**
   * Maximum flavour of incoming quarks
   */
  int maxFlavour_;

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement me_;

  /**
   * Pointer to the H-> 2 gluon vertex (used in gg->H)
   */
  AbstractVVSVertexPtr HGGVertex_;

  /**
   * Pointer to the fermion-fermion Higgs vertex (used in qq->H)
   */
  AbstractFFSVertexPtr HFFVertex_;

  /**
   *  The mass generator for the Higgs
   */
  GenericMassGeneratorPtr hmass_;

  /**
   *  On-shell mass for the higgs
   */
  Energy mh_;

  /**
   *  On-shell width for the higgs
   */
  Energy wh_;

  /**
   *  Stuff for the ME correction
   */
  //@{
  /**
   *  Parameters for the evaluation of the loops for the 
   *  matrix elements
   */
  //@{
  /**
   *  Minimum flavour of quarks to include in the loops
   */
  unsigned int minLoop_;

  /**
   *  Maximum flavour of quarks to include in the loops
   */
  unsigned int maxLoop_;

  /**
   *  Option for treatment of the fermion loops
   */
  unsigned int massOption_;

  /**
   *  Option for dynamic scale choice in alpha_S (0=mT,>0=pT)
   */
  unsigned int mu_R_opt_;

  /**
   *  Option for dynamic scale choice in PDFs    (0=mT,>0=pT)
   */
  unsigned int mu_F_opt_;
  //@}

  //@}

  /**
   *  Small complex number to regularize some integrals
   */
  static const complex<Energy2> epsi_;

  /**
   *  Storage of the loop functions
   */
  //@{
  /**
   *  B functions
   */
  mutable Complex bi_[5];

  /**
   *  C functions
   */
  mutable complex<InvEnergy2> ci_[8];

  /**
   *  D functions
   */
  mutable complex<InvEnergy4> di_[4];
  //@}

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alpha_;

  /**
   *  Mass squared of Higgs
   */  
  Energy2 mh2_;

  /**
   *  Relative weight of the \f$qg\f$ to the \f$gg\f$  channel
   */
  double channelwgtA_;

  /**
   *  Relative weight for the \f$\bar{q}g\f$ to the \f$gg\f$  channel
   */
  double channelwgtB_;

  /**
   *  Weights for the channels as a vector
   */
  vector<double> channelWeights_;

  /**
   *  Power for the \f$\frac{{\rm d}\hat{s}}{\hat{s}^n}\f$ importance sampling
   *  of the \f$gg\f$ component 
   */
  double ggPow_;

  /**
   *  Power for the \f$\frac{{\rm d}\hat{s}}{\hat{s}^n}\f$ importance sampling
   *  of the \f$qg\f$ and \f$\bar{q}g\f$ components 
   */
  double qgPow_;

  /**
   *  The enhancement factor for initial-state radiation
   */
  double enhance_;
  
  /**
   *  Number of weights greater than 1
   */
  unsigned int nover_;

  /**
   *  Number of attempts
   */
  unsigned int ntry_;

  /**
   *  Number which suceed
   */
  unsigned int ngen_;

  /**
   *  Maximum weight
   */
  double maxwgt_;
  //@}

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double power_;

  /**
   *  The prefactor, \f$c\f$ for the \f$gg\f$ channel
   */
  double pregg_;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double preqg_;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double pregqbar_;

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
   *  Properties of the boson and jets
   */
  //@{
  /**
   *  The rapidity of the Higgs boson
   */
  double yh_;

  /**
   *  The mass of the Higgs boson
   */
  Energy mass_;

  /**
   *  the rapidity of the jet
   */
  double yj_;

  /**
   *  The transverse momentum of the jet
   */
  Energy pt_;

  /**
   *  The outgoing parton
   */
  tcPDPtr out_;
  //@}

  /**
   *  Whether of not to construct the vertex for spin correlations
   */
  bool spinCorrelations_;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes of MEPP2Higgs. */
template <>
struct BaseClassTrait<Herwig::MEPP2Higgs,1> {
  /** Typedef of the first base class of MEPP2Higgs. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2Higgs class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2Higgs>
  : public ClassTraitsBase<Herwig::MEPP2Higgs> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2Higgs"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2Higgs is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2Higgs depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2Higgs_H */
