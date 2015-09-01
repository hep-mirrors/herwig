// -*- C++ -*-
//
// MEee2gZ2qq.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_MEee2gZ2qq_H
#define HERWIG_MEee2gZ2qq_H
//
// This is the declaration of the MEee2gZ2qq class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Shower/Couplings/ShowerAlpha.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2gZ2qq class implements the matrix element
 * for \f$e^+e^-\to Z/\gamma \to q\bar{q}\f$ including spin correlations.
 * The class includes greater control over the type of quark produced than is available
 * in the corresponding matrix element from ThePEG, in addition to spin correlations.
 *
 * @see \ref MEee2gZ2qqInterfaces "The interfaces"
 * defined for MEee2gZ2qq.
 */
class MEee2gZ2qq: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEee2gZ2qq() : minflav_(1), maxflav_(5), massopt_(1),
		 spinCorrelations_(true),
		 pTminQED_(GeV), pTminQCD_(GeV),
		 preFactor_(6.)
  {}

  /**
   *  Members for hard corrections to the emission of QCD radiation 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual POWHEGType hasPOWHEGCorrection() {return FSR;}

  /**
   *  Has an old fashioned ME correction
   */
  virtual bool hasMECorrection() {return true;}

  /**
   *  Initialize the ME correction
   */
  virtual void initializeMECorrection(ShowerTreePtr, double &,
				      double & );

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
				     ShowerParticlePtr parent,
				     Branching br);

  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr,
				      vector<ShowerInteraction::Type>);
  //@}

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
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
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
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

protected:

  /**
   *  Calculate the matrix element for \f$e^+e^-\to q \bar{q}\f$.
   * @param partons The incoming and outgoing particles
   * @param momenta The momenta of the incoming and outgoing particles
   * @param first Whether or not to calculate the spin correlations
   */  
  double loME(const vector<cPDPtr> & partons, 
	      const vector<Lorentz5Momentum> & momenta,
	      bool first) const;

  /**
   * Member to calculate the matrix element
   * @param fin  Spinors for incoming fermion
   * @param ain  Spinors for incoming antifermion
   * @param fout Spinors for outgoing fermion
   * @param aout Spinors for outgong antifermion
   * @param me   Spin summed Matrix element
   * @param cont The continuum piece of the matrix element
   * @param BW   The Z piece of the matrix element
   */
  ProductionMatrixElement HelicityME(vector<SpinorWaveFunction>    & fin,
				     vector<SpinorBarWaveFunction> & ain,
				     vector<SpinorBarWaveFunction> & fout,
				     vector<SpinorWaveFunction>    & aout,
				     double & me,
				     double & cont,
				     double & BW ) const;

  /**
   *  The ratio of the matrix element for one additional jet over the
   * leading order result. In practice
   * \f[\frac{\hat{s}|\overline{\mathcal{M}}|^2_2|D_{\rm emit}|}{4\pi C_F\alpha_S|\overline{\mathcal{M}}|^2_3\left(|D_{\rm emit}|+|D_{\rm spect}|\right)}\f]
   * is returned where \f$\|\overline{\mathcal{M}}|^2\f$ is 
   * the spin and colour summed/averaged matrix element.
   * @param partons The incoming and outgoing particles
   * @param momenta The momenta of the incoming and outgoing particles
   * @param iemitter Whether the quark or antiquark is regardede as the emitter
   * @param inter The type of interaction
   * @param subtract Whether or not to subtract the relevant dipole term
   */
  double meRatio(vector<cPDPtr> partons, 
		 vector<Lorentz5Momentum> momenta,
		 unsigned int iemitter,
		 ShowerInteraction::Type inter,
		 bool subtract =false) const;

  /**
   *  Calculate the matrix element for \f$e^-e^-\to q \bar q g\f$.
   * @param partons The incoming and outgoing particles
   * @param momenta The momenta of the incoming and outgoing particles
   * @param inter The type of interaction
   */ 
  InvEnergy2 realME(const vector<cPDPtr> & partons, 
		    const vector<Lorentz5Momentum> & momenta,
		    ShowerInteraction::Type inter) const;

private:

  /**
   *  Generate the momenta for a hard configuration
   */
  pair<Energy,ShowerInteraction::Type> 
  generateHard(ShowerTreePtr tree, 
	       vector<Lorentz5Momentum> & emission,
	       unsigned int & iemit, unsigned int & ispect,
	       bool applyVeto,
	       vector<ShowerInteraction::Type>);

  /**
   *  Calculate \f$\tilde{\kappa}\f$.
   */
  double getKfromX(double, double);

  /**
   *  Vector and axial vector parts of the matrix element
   */
  //@{
  /**
   *  Vector part of the matrix element
   */
  double MEV(double, double);

  /**
   * The matrix element, given \f$x_1\f$, \f$x_2\f$.
   * @param x1 \f$x_1\f$
   * @param x2 \f$x_2\f$
   */
  double PS(double x1, double x2);
  //@}

protected:
  
  /**
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr FFZVertex() const {return FFZVertex_;}
  
  /**
   *  Pointer to the fermion-antifermion photon vertex
   */
  AbstractFFVVertexPtr FFPVertex() const {return FFPVertex_;}
  
  /**
   *  Pointer to the particle data object for the Z
   */
  PDPtr Z0() const {return Z0_;}

  /**
   *  Pointer to the particle data object for the photon
   */
  PDPtr gamma() const {return gamma_;}

  /**
   *  Pointer to the particle data object for the gluon
   */
  PDPtr gluon() const {return gluon_;}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2gZ2qq & operator=(const MEee2gZ2qq &);

private:

  /**
   *  Parameters controlling the leading-order process
   */
  //@{
  /**
   *  The minimum PDG of the quarks to be produced
   */
   int minflav_;

  /**
   *  The maximum PDG of the quarks to be produced
   */
   int maxflav_;

  /**
   *  Option for the treatment of the top quark mass
   */
  unsigned int massopt_;
  //@}

  /**
   *  Pointers to the vertices
   */
  //@{
  /**
   *  Pointer to the fermion-antifermion Z vertex
   */
  AbstractFFVVertexPtr FFZVertex_;
  
  /**
   *  Pointer to the fermion-antifermion photon vertex
   */
  AbstractFFVVertexPtr FFPVertex_;
  
  /**
   *  Pointer to the fermion-antifermion photon vertex
   */
  AbstractFFVVertexPtr FFGVertex_;
  //@}

  /**
   *  Switch on/off the helivity vertex construction
   */
  bool spinCorrelations_;

  /**
   *  Pointer to the ParticleData objects
   */
  //@{
  /**
   *  Pointer to the particle data object for the Z
   */
  PDPtr Z0_;

  /**
   *  Pointer to the particle data object for the photon
   */
  PDPtr gamma_;

  /**
   *  Pointer to the particle data object for the gluon
   */
  PDPtr gluon_;
  //@}

  /**
   * CM energy 
   */
  Energy d_Q_;

  /**
   *  Quark mass
   */
  Energy d_m_;

  /**
   * The rho parameter 
   */
  double d_rho_;

  /**
   * The v parameter
   */
  double d_v_;

  /**
   * The initial kappa-tilde values for radiation from the quark
   */
  double d_kt1_;

  /**
   * The initial kappa-tilde values for radiation from the antiquark
   */
  double d_kt2_;

  /**
   *  Cut-off parameter
   */
  static const double EPS_;

  /**
   *  Pointer to the strong coupling
   */
  ShowerAlphaPtr alphaQCD_;

  /**
   *  Pointer to the EM coupling
   */
  ShowerAlphaPtr alphaQED_;

private:

  /**
   *  Variables for the POWHEG style corrections
   */
  //@{
  /**
   *  The cut off on pt for QED, assuming massless quarks.
   */
  Energy pTminQED_;
  /**
   *  The cut off on pt for QCD, assuming massless quarks.
   */
  Energy pTminQCD_;

  /**
   *  Overestimate for the prefactor
   */
  double preFactor_;

  /**
   *  ParticleData objects for the partons
   */
  vector<cPDPtr> partons_;

  /**
   *  Momenta of the leading-order partons
   */
  vector<Lorentz5Momentum> loMomenta_;
  //@}

};

}

#endif /* HERWIG_MEee2gZ2qq_H */
