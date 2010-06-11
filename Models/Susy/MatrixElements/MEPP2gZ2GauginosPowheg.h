// -*- C++ -*-
#ifndef HERWIG_MEPP2gZ2GauginosPowheg_H
#define HERWIG_MEPP2gZ2GauginosPowheg_H
//
// This is the declaration of the MEPP2gZ2GauginosPowheg class.
//

#include "NLODrellYanBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "Herwig++/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2gZ2GauginosPowheg class implements the matrix elements
 * together with the next-to-leading order corrections for
 * \f$q\bar q\to \gamma Z^0\to \tilde{\ell}\tilde{\ell}^*\f$.
 *
 * @see \ref MEPP2gZ2GauginosPowhegInterfaces "The interfaces"
 * defined for MEPP2gZ2GauginosPowheg.
 */
class MEPP2gZ2GauginosPowheg: public NLODrellYanBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEPP2gZ2GauginosPowheg();
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

  /**
   * Matrix element for \f$q\bar{q}\to \gamma/Z \to f\bar{f}\f$.
   * @param fin  Spinors for incoming quark
   * @param ain  Spinors for incoming antiquark
   * @param fout Spinors for incoming quark
   * @param aout Spinors for incoming antiquark
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double qqbarME(vector<SpinorWaveFunction>    & fin ,
		 vector<SpinorBarWaveFunction> & ain ,
		 vector<SpinorWaveFunction>    & fout ,
		 vector<SpinorBarWaveFunction> & aout ,
		 //ScalarWaveFunction & s1,ScalarWaveFunction &s2,
		 bool me) const;

protected:

  /**
   *  Implementation of the virtual functions for the NLO calculation
   */
  //@{
  /**
   * Virtual matrix element, to be implemented in the
   * inheriting classes. The method should return the
   * loop matrix element including it'sa singular terms
   * Assumes the matrix element has the form
   * \f[ \frac{\alpha_S}{2\pi}\frac{C_F}{\Gamma(1-\epsilon)}
   *     \left(\frac{4\pi\mu^2}{\hat s}\right)^epsilon
   *    \left(\frac{A}{\epsilon^2}+\frac{B}{\epsilon}+C}
   * \f]
   * and A, B and C are returned.
   */


  virtual NLODrellYanBase::Singular virtualME() const;

  /**
   * The leading-order matrix element, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  virtual double loME(const cPDVector & particles,
		      const vector<Lorentz5Momentum> & momenta,
		      bool first=false) const;

  /**
   * The real matrix element divided by \f$2 g_S^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */

  virtual double realME(const cPDVector & particles,
     const vector<Lorentz5Momentum> & momenta) const;
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
  static ClassDescription<MEPP2gZ2GauginosPowheg> initMEPP2gZ2GauginosPowheg;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2gZ2GauginosPowheg & operator=(const MEPP2gZ2GauginosPowheg &);

private:

  /**
   *  Pointer to the vertices for the helicity calculations
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

  /**
   *  Pointer to the chargino Z vertex
   */
  AbstractFFVVertexPtr CCZVertex_;

  /**
   *  Pointer to the neutralino Z vertex
   */
  AbstractFFVVertexPtr NNZVertex_;

  /**
   *  Pointer to the neutralino photon vertex
   */
  AbstractFFVVertexPtr NNPVertex_;

  //@}

  /**
   *  Pointers to the intermediate resonances
   */
  //@{
  /**
   *  Pointer to the Z ParticleData object
   */
  tcPDPtr Z0_;

  /**
   *  Pointer to the photon ParticleData object
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
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2gZ2GauginosPowheg. */
template <>
struct BaseClassTrait<Herwig::MEPP2gZ2GauginosPowheg,1> {
  /** Typedef of the first base class of MEPP2gZ2GauginosPowheg. */
  typedef Herwig::NLODrellYanBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2gZ2GauginosPowheg class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2gZ2GauginosPowheg>
  : public ClassTraitsBase<Herwig::MEPP2gZ2GauginosPowheg> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2gZ2GauginosPowheg"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2gZ2GauginosPowheg is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2gZ2GauginosPowheg depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwSusy.so HwSusyNLO.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2gZ2GauginosPowheg_H */
