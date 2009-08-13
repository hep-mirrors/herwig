// -*- C++ -*-
#ifndef HERWIG_GGtoHHardGenerator_H
#define HERWIG_GGtoHHardGenerator_H
//
// This is the declaration of the GGtoHHardGenerator class.
//

#include "Herwig++/Shower/Base/HardestEmissionGenerator.h"
#include "Herwig++/Shower/Couplings/ShowerAlpha.h"
#include "Herwig++/Shower/Base/ShowerProgenitor.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The GGtoHHardGenerator class implements the generation of hard QCD radiation in
 * \f$gg\to h^0\f$ processes in the POWHEG scheme.
 *
 * @see \ref GGtoHHardGeneratorInterfaces "The interfaces"
 * defined for GGtoHHardGenerator.
 */
class GGtoHHardGenerator: public HardestEmissionGenerator {

  /**
   * Typedef for the BeamParticleData object
   */
  typedef Ptr<BeamParticleData>::transient_const_pointer tcBeamPtr;

public:

  /**
   * The default constructor.
   */
  GGtoHHardGenerator();

  /**
   *  Implementation of virtual members from HardestEmissionGenerator
   */
  //@{
  /**
   *  Member to generate the hardest emission
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);

  /**
   *  Member to decide if the inheriting class can handle this process
   */
  virtual bool canHandle(ShowerTreePtr);
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
  inline virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:
 
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

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GGtoHHardGenerator> initGGtoHHardGenerator;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GGtoHHardGenerator & operator=(const GGtoHHardGenerator &);

private:

  /**
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr _alphaS;

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   * The power, \f$n\f$, for the sampling
   */
  double _power;

  /**
   *  The prefactor, \f$c\f$ for the \f$gg\f$ channel
   */
  double _pregg;

  /**
   *  The prefactor, \f$c\f$ for the \f$qg\f$ channel
   */
  double _preqg;

  /**
   *  The prefactor, \f$c\f$ for the \f$g\bar{q}\f$ channel
   */
  double _pregqbar;

  /**
   *  The prefactors as a vector for easy use
   */
  vector<double> _prefactor;
  //@}

  /**
   *  The transverse momentum of the jet
   */
  Energy _min_pt;

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

  /**
   *  Prefactor if variable scale used
   */
  double scaleFact_;

  /**
   *  Option for using pt or mt as alpha_S scale
   */
  unsigned int alphaScale_;

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
   *  Properties of the incoming particles
   */
  //@{
  /**
   *  Pointers to the BeamParticleData objects
   */
  vector<tcBeamPtr> _beams;
  
  /**
   *  Pointers to the ParticleDataObjects for the partons
   */
  vector<tcPDPtr> _partons;
  //@}

  /**
   *  Properties of the boson and jets
   */
  //@{
  /**
   *  The rapidity of the Higgs boson
   */
  double _yh;

  /**
   *  The mass of the Higgs boson
   */
  Energy _mass;

  /**
   *  Mass squared of Higgs
   */  
  Energy2 _mh2;

  /**
   *  the rapidity of the jet
   */
  double _yj;

  /**
   *  The transverse momentum of the jet
   */
  Energy _pt;

  /**
   *  The outgoing parton
   */
  tcPDPtr _out;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GGtoHHardGenerator. */
template <>
struct BaseClassTrait<Herwig::GGtoHHardGenerator,1> {
  /** Typedef of the first base class of GGtoHHardGenerator. */
  typedef Herwig::HardestEmissionGenerator NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GGtoHHardGenerator class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::GGtoHHardGenerator>
  : public ClassTraitsBase<Herwig::GGtoHHardGenerator> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::GGtoHHardGenerator"; }
  /**
   * The name of a file containing the dynamic library where the class
   * GGtoHHardGenerator is implemented. It may also include several, space-separated,
   * libraries if the class GGtoHHardGenerator depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegME.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_GGtoHHardGenerator_H */
