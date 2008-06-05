// -*- C++ -*-
//
// MEPP2HiggsJet.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2007 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEPP2HiggsJet_H
#define HERWIG_MEPP2HiggsJet_H
//
// This is the declaration of the MEPP2HiggsJet class.
//

#include "ThePEG/MatrixElement/ME2to2Base.h"
#include "Herwig++/Utilities/Maths.h"
#include "Herwig++/MatrixElement/General/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig++/PDT/SMHiggsMassGenerator.h"
#include "Herwig++/MatrixElement/General/ProductionMatrixElement.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

/**
 * The MEPP2HiggsJet class implements the matrix element for Higgs+jet production.
 *
 * @see \ref MEPP2HiggsJetInterfaces "The interfaces"
 * defined for MEPP2HiggsJet.
 */
class MEPP2HiggsJet: public ME2to2Base {

public:

  /**
   * The default constructor.
   */
  inline MEPP2HiggsJet();

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the matrix element for the kinematical configuation
   * previously provided by the last call to setKinematics(). Uses
   * me().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * The number of internal degreed of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
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
  //@}

private:

  /**
   *  Members to return the matrix elements for the different subprocesses
   */
  //@{
  /**
   * Matrix element for \f$q\bar{q}\to Hg\f$.
   * @param fin   Spinors for incoming quark
   * @param ain   Spinors for incoming antiquark
   * @param hout  Wavefunction for the outgoing higgs
   * @param gout  Polarization vectors for the outgoing gluon
   * @param me    Whether or not to calculate the matrix element for spin correlations
   **/
  double qqbarME(vector<SpinorWaveFunction> & fin, vector<SpinorBarWaveFunction> & ain,
		 ScalarWaveFunction & hout, vector<VectorWaveFunction> & gout,
		 bool me) const;

  /**
   * Matrix element for \f$qg\to Hq\f$.
   * @param fin  Spinors for incoming quark
   * @param gin  Polarization vectors for the incoming gluon
   * @param hout Wavefunction for the outgoing higgs
   * @param fout Spinors for outgoing quark
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
  double qgME(vector<SpinorWaveFunction> & fin,vector<VectorWaveFunction> & gin,
	      ScalarWaveFunction & hout, vector<SpinorBarWaveFunction> & fout,
	      bool me) const;

  /**
   * Matrix element for \f$\bar{q}g\to H\bar{q}\f$.
   * @param fin  Spinors for incoming antiquark
   * @param gin  Polarization vectors for the incoming gluon
   * @param hout Wavefunction for the outgoing higgs
   * @param fout Spinors for outgoing antiquark
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
  double qbargME(vector<SpinorBarWaveFunction> & fin,vector<VectorWaveFunction> & gin,
		 ScalarWaveFunction & hout, vector<SpinorWaveFunction> & fout,
		 bool me) const;

  /**
   * Matrix element for \f$gg\to Hg\f$.
   * @param g1   Polarization vectors for the first  incoming gluon
   * @param g2   Polarization vectors for the second incoming gluon
   * @param hout Wavefunction for the outgoing higgs
   * @param g4   Polarization vectors for the outgoing gluon
   * @param me   Whether or not to calculate the matrix element for spin correlations
   **/
  double ggME(vector<VectorWaveFunction> g1, vector<VectorWaveFunction> g2,
	      ScalarWaveFunction & hout,     vector<VectorWaveFunction> g4,
	      bool me) const;
  //@}

private:

  /**
   *  Members to calculate the functions for the loop diagrams
   */
  //@{
  /**
   *  The \f$W_1(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  inline Complex W1(Energy2 s,Energy2 mf2) const;

  /**
   *  The \f$W_2(s)\f$ function of NPB297 (1988) 221-243.
   * @param s   The invariant
   * @param mf2 The fermion mass squared
   */
  inline Complex W2(Energy2 s,Energy2 mf2) const;

  /**
   * The \f$W_3(s,t,u,v)\f$ function of NPB297 (1988) 221-243.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param v   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  inline Complex W3(Energy2 s, Energy2 t, Energy2 u, Energy2 v, Energy2 mf2) const;

  /**
   * The \f$I_3(s,t,u,v)\f$ function of NPB297 (1988) 221-243.
   * @param s The \f$s\f$ invariant
   * @param t The \f$t\f$ invariant
   * @param u The \f$u\f$ invariant
   * @param v The \f$v\f$ invariant
   * @param mf2 The fermion mass squared
   */
  inline Complex I3(Energy2 s, Energy2 t, Energy2 u, Energy2 v, Energy2 mf2) const;

  /**
   * The \f$b_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  inline Complex b2(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) const;

  /**
   * The \f$b_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  inline Complex b4(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) const;

  /**
   * The \f$A_2(s,t,u)\f$ function of NPB297 (1988) 221-243.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  inline Complex A2(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) const;

  /**
   * The \f$A_4(s,t,u)\f$ function of NPB297 (1988) 221-243.
   * @param s   The \f$s\f$ invariant
   * @param t   The \f$t\f$ invariant
   * @param u   The \f$u\f$ invariant
   * @param mf2 The fermion mass squared.
   */
  inline Complex A4(Energy2 s, Energy2 t, Energy2 u, Energy2 mf2) const;
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2HiggsJet> initMEPP2HiggsJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2HiggsJet & operator=(const MEPP2HiggsJet &);

private:

  /**
   * Defines the Higgs resonance shape
   */
  unsigned int _shapeopt;

  /**
   *  Maximum PDG code of the quarks allowed
   */
  unsigned int _maxflavour;

  /**
   *  Option for which processes to include
   */
  unsigned int _process;
  
  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement _me;

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
   *  Small complex number to regularize some integrals
   */
  static const Complex _epsi;

  /**
   *  On-shell mass for the higgs
   */
  Energy _mh;

  /**
   *  On-shell width for the higgs
   */
  Energy _wh;

  /**
   *  The mass generator for the Higgs
   */
  SMHiggsMassGeneratorPtr _hmass;

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
  mutable Complex _ci[8];

  /**
   *  D functions
   */
  mutable Complex _di[4];
  //@}

  /**
   *  Storage of the diagram weights for the \f$gg\to Hg\f$ subprocess
   */
  mutable double _diagwgt[3];
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2HiggsJet. */
template <>
struct BaseClassTrait<Herwig::MEPP2HiggsJet,1> {
  /** Typedef of the first base class of MEPP2HiggsJet. */
  typedef ME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2HiggsJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2HiggsJet>
  : public ClassTraitsBase<Herwig::MEPP2HiggsJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2HiggsJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2HiggsJet is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2HiggsJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#include "MEPP2HiggsJet.icc"

#endif /* HERWIG_MEPP2HiggsJet_H */
