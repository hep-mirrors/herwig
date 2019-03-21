// -*- C++ -*-
#ifndef HERWIG_MEPP2QQHiggs_H
#define HERWIG_MEPP2QQHiggs_H
//
// This is the declaration of the MEPP2QQHiggs class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/PDT/GenericMassGenerator.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2QQHiggs class implements the matrix elements for
 * \f$gg\to Q \bar Q h^0\f$ and \f$q\bar q\to Q \bar Q h^0\f$.
 *
 * @see \ref MEPP2QQHiggsInterfaces "The interfaces"
 * defined for MEPP2QQHiggs.
 */
class MEPP2QQHiggs: public HwMEBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MEPP2QQHiggs();
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
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

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

  /**
   *  Construct the vertex of spin correlations.
   */
  virtual void constructVertex(tSubProPtr);
  //@}

protected:

  /**
   *  Members to calculate the matrix elements
   */
  //@{
  /**
   * Matrix element for \f$gg\to Q\bar{Q}h^0\f$
   * @param g1   The wavefunctions for the first  incoming gluon
   * @param g2   The wavefunctions for the second incoming gluon
   * @param q    The wavefunction  for the outgoing quark
   * @param qbar The wavefunction  for the outgoing antiquark
   * @param h    The wavefunction for the outgoing Higgs boson
   * @param flow The colour flow
   */
  double ggME(vector<VectorWaveFunction> &g1,vector<VectorWaveFunction> &g2,
	      vector<SpinorBarWaveFunction> & q,vector<SpinorWaveFunction> & qbar,
	      ScalarWaveFunction & h,
	      unsigned int flow) const;

  /**
   * Matrix element for \f$q\bar{q}\to Q\bar{Q}h^0\f$
   * @param q1 The wavefunction  for the incoming quark
   * @param q2 The wavefunction  for the incoming antiquark
   * @param q3 The wavefunction  for the outgoing quark
   * @param q4 The wavefunction  for the outgoing antiquark
   * @param h    The wavefunction for the outgoing Higgs boson
   * @param flow The colour flow
   */
  double qqME(vector<SpinorWaveFunction> & q1,
	      vector<SpinorBarWaveFunction> & q2,
	      vector<SpinorBarWaveFunction>    & q3,
	      vector<SpinorWaveFunction>    & q4,
	      ScalarWaveFunction & h,
	      unsigned int flow) const;
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
  static ClassDescription<MEPP2QQHiggs> initMEPP2QQHiggs;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2QQHiggs & operator=(const MEPP2QQHiggs &) = delete;

private:
  
  /**
   *  Switches to control the subprocess
   */
  //@{
  /**
   *  Quark Flavour
   */
  unsigned int quarkFlavour_;
  
  /**
   *  Processes to include
   */
  unsigned int process_;
  //@}

  /**
   *  Switches etc for the Higgs mass generation
   */
  //@{
  /**
   * Defines the Higgs resonance shape
   */
  unsigned int shapeOpt_;

  /**
   *  On-shell mass for the higgs
   */
  Energy mh_;

  /**
   *  On-shell width for the higgs
   */
  Energy wh_;

  /**
   *  The mass generator for the Higgs
   */
  GenericMassGeneratorPtr hmass_;
  //@}
  
  /**
   *  Vertices needed to compute the diagrams
   */
  //@{
  /**
   *  \f$ggg\f$ vertex
   */
  AbstractVVVVertexPtr GGGVertex_;
  
  /**
   *  \f$q\bar{q}g\f$ vertex
   */
  AbstractFFVVertexPtr QQGVertex_;

  /**
   *  \f$q\bar q h^0\f$ vertex
   */
  AbstractFFSVertexPtr QQHVertex_;
  //@}
  

  /**
   *  ParticleData objects of the particles
   */  
  //@{
  /**
   *  The gluon
   */
  PDPtr gluon_;

  /**
   * The Higgs boson
   */
  PDPtr higgs_;
  
  /**
   *  the quarks
   */
  vector<PDPtr> quark_;

  /**
   *  the antiquarks
   */
  vector<PDPtr> antiquark_;
  //@}

  /**
   *  Parameters for the phase-space generation
   */
  //@{
  /**
   *  Power for the phase-space mapping
   */
  double alpha_;
  //@}

  /**
   *  Colour flow
   */
  mutable unsigned int flow_;

  /**
   *  Diagram
   */
  mutable unsigned int diagram_;

  /**
   *  Matrix element
   */
  mutable ProductionMatrixElement me_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2QQHiggs. */
template <>
struct BaseClassTrait<Herwig::MEPP2QQHiggs,1> {
  /** Typedef of the first base class of MEPP2QQHiggs. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2QQHiggs class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2QQHiggs>
  : public ClassTraitsBase<Herwig::MEPP2QQHiggs> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2QQHiggs"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2QQHiggs is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2QQHiggs depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2QQHiggs_H */
