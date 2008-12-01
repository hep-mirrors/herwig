// -*- C++ -*-
#ifndef HERWIG_GGtoHMECorrection_H
#define HERWIG_GGtoHMECorrection_H
//
// This is the declaration of the GGtoHMECorrection class.
//

#include "QTildeMECorrection.h"
#include "Herwig++/Utilities/Maths.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GGtoHMECorrection class implements the matrix element correction
 * for \f$gg\to h\f$.
 *
 * @see \ref GGtoHMECorrectionInterfaces "The interfaces"
 * defined for GGtoHMECorrection.
 */
class GGtoHMECorrection: public QTildeMECorrection {

/**
 *  Typedef for BeamParticleData pointers
 */
typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  GGtoHMECorrection();

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
  bool applyHard(ShowerParticleVector gluons,
		 vector<tcBeamPtr> beams,
		 PPtr higgs,unsigned int & iemit,
		 unsigned int & itype,vector<Lorentz5Momentum> & pnew,
		 pair<double,double> & xnew,
		 tPDPtr & out);

  /**
   *   Members to calculate the matrix elements
   */
  //@{
  /**
   *  The leading-order matrix element for \f$gg\to H\f$
   */
  Energy4 loME();

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
  Complex F(double x);
  //@}

  /**
   *  Method to extract the PDF weight for quark/antiquark
   * initiated processes and select the quark flavour
   */  
  tPDPtr quarkFlavour(tcPDFPtr pdf, Energy2 scale, double x, tcBeamPtr beam, 
		      double & pdfweight, bool anti);

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
  static ClassDescription<GGtoHMECorrection> initGGtoHMECorrection;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GGtoHMECorrection & operator=(const GGtoHMECorrection &);

private:

  /**
   *  Parameters for the evaluation of the loops for the 
   *  matrix elements
   */
  //@{
  /**
   *  Minimum flavour of quarks to include in the loops
   */
  unsigned int _minloop;

  /**
   *  Maximum flavour of quarks to include in the loops
   */
  unsigned int _maxloop;

  /**
   *  Option for treatment of the fermion loops
   */
  unsigned int _massopt;
  //@}

  /**
   *  Small complex number to regularize some integrals
   */
  static const complex<Energy2> _epsi;

  /**
   *  Storage of the loop functions
   */
  //@{
  /**
   *  B functions
   */
  mutable Complex _bi[5];

  /**
   *  C functions
   */
  mutable complex<InvEnergy2> _ci[8];

  /**
   *  D functions
   */
  mutable complex<InvEnergy4> _di[4];
  //@}

  /**
   *  Storage of the diagram weights for the \f$gg\to Hg\f$ subprocess
   */
  mutable double _diagwgt[3];

  /**
   *  Mass squared of Higgs
   */  
  Energy2 _mh2;

  /**
   *  Relative weight of the \f$qg\f$ to the \f$gg\f$  channel
   */
  double _channelwgtA;

  /**
   *  Relative weight for the \f$\bar{q}g\f$ to the \f$gg\f$  channel
   */
  double _channelwgtB;

  /**
   *  Weights for the channels as a vector
   */
  vector<double> _channelweights;

  /**
   *  Power for the \f$\frac{{\rm d}\hat{s}}{\hat{s}^n}\f$ importance sampling
   *  of the \f$gg\f$ component 
   */
  double _ggpow;

  /**
   *  Power for the \f$\frac{{\rm d}\hat{s}}{\hat{s}^n}\f$ importance sampling
   *  of the \f$qg\f$ and \f$\bar{q}g\f$ components 
   */
  double _qgpow;

  /**
   *  The enhancement factor for initial-state radiation
   */
  double _enhance;
  
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
  double _maxwgt;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GGtoHMECorrection. */
template <>
struct BaseClassTrait<Herwig::GGtoHMECorrection,1> {
  /** Typedef of the first base class of GGtoHMECorrection. */
  typedef Herwig::QTildeMECorrection NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GGtoHMECorrection class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GGtoHMECorrection>
  : public ClassTraitsBase<Herwig::GGtoHMECorrection> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GGtoHMECorrection"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GGtoHMECorrection is implemented. It may also include several, space-separated,
   * libraries if the class GGtoHMECorrection depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_GGtoHMECorrection_H */
