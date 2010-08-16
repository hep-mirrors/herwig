// -*- C++ -*-
#ifndef HERWIG_MEqq2gZ2ffPowhegQED_H
#define HERWIG_MEqq2gZ2ffPowhegQED_H
//
// This is the declaration of the MEqq2gZ2ffPowhegQED class.
//

#include "Herwig++/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEqq2gZ2ffPowhegQED class implements 
 * \f$q\bar q \to \gamma/Z^0 \to \ell^+\ell^-\f$
 * including both QED and QCD corrections
 *
 * @see \ref MEqq2gZ2ffPowhegQEDInterfaces "The interfaces"
 * defined for MEqq2gZ2ffPowhegQED.
 */
class MEqq2gZ2ffPowhegQED: public HwMEBase {

public:

  /**
   *  Enum for the type of Dipole
   */
  enum DipoleType {
    II12,II21,FF34,FF43
  };

public:

  /**
   * The default constructor.
   */
  MEqq2gZ2ffPowhegQED();

  /**
   *  Virtual members to be overridden by inheriting classes
   *  which implement hard corrections 
   */
  //@{
  /**
   *  Has a POWHEG style correction
   */
  virtual bool hasPOWHEGCorrection() {return true;}
  
  /**
   *  Apply the POWHEG style correction
   */
  virtual HardTreePtr generateHardest(ShowerTreePtr);
  //@}

public:

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

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;
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
   *  Virtual functions for the NLO calculation
   */
  //@{
  /**
   *  Calculate of the full next-to-leading order weight
   */
  double NLOWeight() const;

  /**
   * The leading-order matrix element, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double loME(const cPDVector & particles,
	      const vector<Lorentz5Momentum> & momenta,
	      bool first=false) const;
  
  /**
   * The real matrix element divided by \f$2 g_S^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  double realQCDME(const cPDVector & particles,
		   const vector<Lorentz5Momentum> & momenta) const;
  
  /**
   * The real matrix element divided by \f$2 e^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  vector<double> realQEDME(const cPDVector & particles,
			   const vector<Lorentz5Momentum> & momenta) const;
  //@}

  /**
   * Matrix element for \f$q\bar{q}\to \gamma/Z \to f\bar{f}\f$.
   * @param fin  Spinors for incoming quark
   * @param ain  Spinors for incoming antiquark
   * @param fout Spinors for outgoing lepton
   * @param aout Spinors for outgoing antilepton
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction>    & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<SpinorBarWaveFunction> & fout,
		 vector<SpinorWaveFunction>    & aout,
		 bool me) const;


  /**
   *  Subtracted virtual contribution
   */
  double subtractedVirtual() const;

  /**
   *  Subtracted real contribution
   */
  vector<double> 
  subtractedRealQCD(pair<double,double> x, double z,
		    double zJac,
		    double oldqPDF, double newqPDF, double newgPDF,
		    DipoleType dipole) const;

  /**
   *  Subtracted real contribution
   */
  vector<double> 
  subtractedRealQED(pair<double,double> x, double z,
		    double zJac,
		    double oldqPDF, double newqPDF, double newpPDF,
		    DipoleType dipole) const;

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
   *  Gluon or photon collinear counter term
   */
  double collinearBoson(Energy2 mu2, double jac, double z,
			double oldPDF, double newPDF) const;
  //@}

  /**
   *   Dipole subtracted matrix elements
   */
  //@{
  /**
   *  \f$q\bar q\f$
   */
  pair<double,double> 
  subtractedQCDMEqqbar(const vector<Lorentz5Momentum> & pnew, 
		       DipoleType dipole, bool subtract) const;
  
  /**
   *  \f$g\bar q\f$
   */
  pair<double,double> 
  subtractedQCDMEgqbar(const vector<Lorentz5Momentum> & pnew,
		       DipoleType dipole, bool subtract) const;
  
  /**
   *  \f$q\bar q\f$
   */
  pair<double,double> 
  subtractedQEDMEqqbar(const vector<Lorentz5Momentum> & pnew,
		       DipoleType dipole, bool subtract) const;
  
  /**
   *  \f$g\bar q\f$
   */
  pair<double,double> 
  subtractedQEDMEpqbar(const vector<Lorentz5Momentum> & pnew,
		       DipoleType dipole, bool subtract) const;

  /**
   *  The supression function
   */
  pair<double,double> supressionFunction(Energy2 pt2) const {
    switch( supressionFunction_ ) {
    case 0:
      return make_pair(1.,0.);
    case 1:
      if(pt2 <  lambda2_ ) return make_pair(1.,0.);
      else                 return make_pair(0.,1.);
    case 2:
      return make_pair(lambda2_/(pt2+lambda2_),pt2/(pt2+lambda2_));
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
   * The static object used to initialize the description of this class.
   * Indicates that this is a class with persistent data.
   */
  static ClassDescription<MEqq2gZ2ffPowhegQED> initMEqq2gZ2ffPowhegQED;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEqq2gZ2ffPowhegQED & operator=(const MEqq2gZ2ffPowhegQED &);

private:

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
   *   Which corrections to included
   */
  unsigned int corrections_;

  /**
   *  Include incoming photons?
   */
  bool incomingPhotons_;

  /**
   *  Whether to generate the positive, negative or leading order contribution
   */
  unsigned int contrib_;

  /**
   *  Power for sampling \f$x_p\f$
   */
  double power_;

  /**
   *  Phase-space sampling for z
   */
  double zPow_;

  /**
   *  Phase-space sampling for y
   */
  double yPow_;

  /**
   *  Pointer to the gluon ParticleData object
   */
  tcPDPtr gluon_;

  /**
   *  Factor for \f$C_F\f$ dependent pieces
   */
  mutable double CFfact_;

  /**
   *  Factor for \f$T_R\f$ dependent pieces
   */
  mutable double TRfact_;

  /**
   *  Factor for EM pieces
   */
  mutable double EMfact_;

  /**
   *  Strong coupling
   */
  mutable double alphaS_;

  /**
   *  Strong coupling
   */
  mutable double alphaEM_;

  /**
   *  Use a fixed value of couplings
   */
  bool fixedCouplings_;

  /**
   *  Leading-order matrix element
   */
  mutable double loME_;

  /**
   *  Structure for the virtual FS QED corrections
   */
  mutable double f2term_;

  /**
   *  Choice of the supression function
   */
  unsigned int supressionFunction_;

  /**
   *  Scalar for the supression function
   */
  Energy2 lambda2_;

  /**
   *  Storage of the weights of the different processes for the hard emission
   * generation
   */
  mutable vector<double> weights_;

private:

  /**
   *  Momenta of the particles for gluon emmision from first the particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQCDGluon1_;

  /**
   *  Momenta of the particles for anti quark emission from the first particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQCDQuark1_;

  /**
   *  Momenta of the particle for gluon emission from the second particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQCDGluon2_;

  /**
   *  Momenta of the particles for quark emission from the second particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQCDQuark2_;

  /**
   *  Momenta of the particles for gluon emmision from first the particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQEDPhoton1_;

  /**
   *  Momenta of the particles for anti quark emission from the first particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQEDQuark1_;

  /**
   *  Momenta of the particle for gluon emission from the second particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQEDPhoton2_;

  /**
   *  Momenta of the particles for quark emission from the second particle
   */
  mutable vector<Lorentz5Momentum> realEmissionQEDQuark2_;

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
   *  Whether the quark is in the + or - z direction
   */
  bool _quarkplus;

  /**
   *  Constants for the sampling. The distribution is assumed to have the
   *  form \f$\frac{c}{{\rm GeV}}\times\left(\frac{{\rm GeV}}{p_T}\right)^n\f$ 
   */
  //@{
  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preqqbarq_;
  /**
   *  The prefactor, \f$c\f$ for the \f$q\bar{q}\f$ channel
   */
  double preqqbarqbar_;

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
   *  Pointer to the object calculating the strong coupling
   */
  ShowerAlphaPtr alphaQCD_;

  /**
   *  Pointers to the intermediate resonances
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr Z0_;

  /**
   *  Pointer to the photon PartcileData object
   */
  tcPDPtr gamma_;
  //@}

  /**
   * Matrix element for spin correlations
   */
  ProductionMatrixElement me_;

  /**
   *  Process
   */
  int process_;

  /**
   *  Maximum flavour of the incoming quarks
   */
  int maxFlavour_;

  /**
   *  Pointers to the vertices for the helicity calculations
   */
  //@{
  /**
   *  Pointer to the Z fermions vertex
   */
  AbstractFFVVertexPtr FFZVertex_;
  
  /**
   *  Pointer to the photon fermions vertex
   */
  AbstractFFVVertexPtr FFPVertex_;
  
  /**
   *  Pointer to the gluon fermions vertex
   */
  AbstractFFVVertexPtr FFGVertex_;
  //@}

};
}


#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEqq2gZ2ffPowhegQED. */
template <>
struct BaseClassTrait<Herwig::MEqq2gZ2ffPowhegQED,1> {
  /** Typedef of the first base class of MEqq2gZ2ffPowhegQED. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEqq2gZ2ffPowhegQED class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEqq2gZ2ffPowhegQED>
  : public ClassTraitsBase<Herwig::MEqq2gZ2ffPowhegQED> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEqq2gZ2ffPowhegQED"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEqq2gZ2ffPowhegQED is implemented. It may also include several, space-separated,
   * libraries if the class MEqq2gZ2ffPowhegQED depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEqq2gZ2ffPowhegQED_H */


