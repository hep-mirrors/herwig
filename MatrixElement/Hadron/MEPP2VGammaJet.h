// -*- C++ -*-
#ifndef HERWIG_MEPP2VGammaJet_H
#define HERWIG_MEPP2VGammaJet_H
//
// This is the declaration of the MEPP2VGammaJet class.
//

#include "Herwig++/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"

namespace Herwig {

using namespace ThePEG;

/**
 * Here is the documentation of the MEPP2VGammaJet class.
 *
 * @see \ref MEPP2VGammaJetInterfaces "The interfaces"
 * defined for MEPP2VGammaJet.
 */
class MEPP2VGammaJet: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2VGammaJet();

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
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /**
   * The real matrix element divided by \f$2 g_S^2\f$, to be implemented in the
   * inheriting classes. 
   * @param particles The ParticleData objects of the particles
   * @param momenta The momenta of the particles
   */
  double realME(const cPDVector & particles,
		const vector<Lorentz5Momentum> & momenta) const;
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
  static ClassDescription<MEPP2VGammaJet> initMEPP2VGammaJet;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VGammaJet & operator=(const MEPP2VGammaJet &);

private:

  /**
   *  The allowed processes
   */
  unsigned int process_;

  /**
   *  The maximum flavour
   */
  unsigned int maxFlavour_;

  /**
   *  Boson
   */
  unsigned int boson_;

  /**
   *  Choice of the phase-space mapping
   */
  unsigned int phase_;

  /**
   *  Strong coupling
   */
  mutable double alphaS_;

  /**
   *  Use a fixed value of \f$\alpha_S\f$
   */
  bool fixedAlphaS_;

  /**
   *  Electromagnetic coupling
   */
  mutable double alphaEM_;

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
  FFVVertexPtr FFWvertex_;

  /**
   *   FFGVertex
   */
  AbstractFFVVertexPtr FFGvertex_;

  /**
   *   FFZVertex
   */
  FFVVertexPtr FFZvertex_;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr WWWvertex_;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VGammaJet. */
template <>
struct BaseClassTrait<Herwig::MEPP2VGammaJet,1> {
  /** Typedef of the first base class of MEPP2VGammaJet. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VGammaJet class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VGammaJet>
  : public ClassTraitsBase<Herwig::MEPP2VGammaJet> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VGammaJet"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VGammaJet is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VGammaJet depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VGammaJet_H */
