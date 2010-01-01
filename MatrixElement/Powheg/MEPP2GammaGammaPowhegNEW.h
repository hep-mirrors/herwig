// -*- C++ -*-
#ifndef HERWIG_MEPP2GammaGammaPowhegNEW_H
#define HERWIG_MEPP2GammaGammaPowhegNEW_H
//
// This is the declaration of the MEPP2GammaGammaPowhegNEW class.
//

#include "Herwig++/MatrixElement/HwME2to2Base.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"
#include "Herwig++/Shower/Base/HardTree.fh"
#include "Herwig++/Shower/Base/ShowerTree.fh"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2GammaGammaPowhegNEW class.
 *
 * @see \ref MEPP2GammaGammaPowhegNEWInterfaces "The interfaces"
 * defined for MEPP2GammaGammaPowhegNEW.
 */
class MEPP2GammaGammaPowhegNEW: public HwME2to2Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEPP2GammaGammaPowhegNEW();
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

  /**
   *  Apply the POWHEG style correction
   */
  HardTreePtr generateHardest(ShowerTreePtr) const;

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
   *  Leading-order matrix elements for photon pair production
   */
  //@{
  /**
   * The leading-order matrix element, for \f$q\bar q \to \gamma \gamma\f$
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double qqbarGammaGammaME(const cPDVector & particles,
			   const vector<Lorentz5Momentum> & momenta,
			   bool first=false) const;
  
  /**
   * The leading-order matrix element, for \f$g g \to \gamma \gamma\f$
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double ggGammaGammaME(const cPDVector & particles,
			const vector<Lorentz5Momentum> & momenta,
			bool first=false) const;

  /**
   *  \f$gg\to\gamma\gamma\f$ matrix element for the \f$++++\f$ helicity configuration.
   * @param s The \f$s\f$ invariant
   * @param t The \f$t\f$ invariant
   * @param u The \f$u\f$ invariant
   */
  Complex ggme(Energy2 s,Energy2 t,Energy2 u) const {
    double ltu(log(abs(t/u)));
    double frac1((t-u)/s),frac2((sqr(t)+sqr(u))/sqr(s));
    double thetatu = (t/u<0) ? 0 : 1;
    double thetat  = (t<ZERO)   ? 0 : 1;
    double thetau  = (u<ZERO)   ? 0 : 1;
    using Constants::pi;
    return Complex(1.+frac1*ltu+0.5*frac2*(sqr(ltu)+sqr(pi)*thetatu),
		   -pi*(thetat-thetau)*(frac1+frac2*ltu));
  }
  //@}

  /**
   *  Leading-order matrix elements for photon+jet production
   */
  //@{
  /**
   * The leading-order matrix element, for \f$q\bar q \to \gamma g\f$
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double qqbarGammaJetME(const cPDVector & particles,
			 const vector<Lorentz5Momentum> & momenta,
			 bool first=false) const;

  /**
   * The leading-order matrix element, for \f$qg  \to \gamma q\f$
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double qgGammaJetME(const cPDVector & particles,
			 const vector<Lorentz5Momentum> & momenta,
			 bool first=false) const;

  /**
   * The leading-order matrix element, for \f$g \bar q \to \gamma q\f$
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   * @param first Whether or not this is the first call and the spin
   * and diagram information should be stored
   */
  double gqbarGammaJetME(const cPDVector & particles,
			 const vector<Lorentz5Momentum> & momenta,
			 bool first=false) const;
  //@}

  /**
   *  Real Emission matrix elements
   */
  //@{
  /**
   *  Matrix element for \f$q\bar q \to \gamma\gamma g\f$
   */
  InvEnergy2 qqbargME(const vector<tcPDPtr> & particles,
		      const vector<Lorentz5Momentum> & p) const;

  /**
   *  Matrix element for \f$q g \to \gamma\gamma q\f$
   */
  InvEnergy2 qgqME(const vector<tcPDPtr> & particles,
		   const vector<Lorentz5Momentum> & p) const;

  /**
   *  Matrix element for \f$g \bar q\to \gamma\gamma \bar q\f$
   */
  InvEnergy2 qbargqbarME(const vector<tcPDPtr> & particles,
			 const vector<Lorentz5Momentum> & p) const;
  //@}

  /**
   *  Supression factors
   */
  //@{
  
  //@}

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
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MEPP2GammaGammaPowhegNEW> initMEPP2GammaGammaPowhegNEW;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2GammaGammaPowhegNEW & operator=(const MEPP2GammaGammaPowhegNEW &);

private:

  /**
   *
   */
  //@{
  /**
   *  Pointer to the quark-antiquark-gluon vertex
   */
  AbstractFFVVertexPtr gluonVertex_;

  /**
   *  Pointer to the quark-antiquark-photon vertex
   */
  AbstractFFVVertexPtr photonVertex_;
  //@}

  /**
   *  Switches controlling which processes to generate
   */
  //@{
  /**
   *  Which processes/order to include
   */
  unsigned int process_;

  /**
   *  Maximum flavour of quarks to include
   */
  unsigned int maxFlavour_;
  //@}

  /**
   *  Parameters controlling the \f$F\f$ functions supressing high \f$p_T\f$
   *  radiation
   */
  //@{
  /**
   *  Option for the QCD supression function
   */
  unsigned int QCDOption_;

  /**
   *  Fixed scale for the QCD factor
   */
  Energy2 QCDScale_;

  /**
   *  Option for the QED supression function
   */
  unsigned int QEDOption_;

  /**
   *  Fixed scale for the QED factor
   */
  Energy2 QEDScale_;
  //@}

  /**
   *  Radiative variables
   */
  //@{
  /**
   *  \f$\tilde{z}\f$
   */
  double ztilde_;

  /**
   *   \f$\tilde{v}\f$
   */
  double vtilde_;

  /**
   *  \f$\phi\f$
   */
  double phi_;
  //@}
  
  /**
   * Matrix element for spin correlations
   */
  mutable ProductionMatrixElement _me;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2GammaGammaPowhegNEW. */
template <>
struct BaseClassTrait<Herwig::MEPP2GammaGammaPowhegNEW,1> {
  /** Typedef of the first base class of MEPP2GammaGammaPowhegNEW. */
  typedef Herwig::HwME2to2Base NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2GammaGammaPowhegNEW class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2GammaGammaPowhegNEW>
  : public ClassTraitsBase<Herwig::MEPP2GammaGammaPowhegNEW> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2GammaGammaPowhegNEW"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2GammaGammaPowhegNEW is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2GammaGammaPowhegNEW depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so HwPowhegMEHadron.so HwPowhegShower.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2GammaGammaPowhegNEW_H */
