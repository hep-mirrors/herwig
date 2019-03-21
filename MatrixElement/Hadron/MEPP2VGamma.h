// -*- C++ -*-
#ifndef HERWIG_MEPP2VGamma_H
#define HERWIG_MEPP2VGamma_H
//
// This is the declaration of the MEPP2VGamma class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEPP2VGamma class implements the .
 *
 * @see \ref MEPP2VGammaInterfaces "The interfaces"
 * defined for MEPP2VGamma.
 */
class MEPP2VGamma: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEPP2VGamma();

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
  
  /**
   * Matrix element for \f$f\bar{f}\to W^\pm \gamma\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  \f$W^\pm\f$ wavefunction
   * @param v2  \f$\gamma\f$ wavefunction
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double WGammaME(vector<SpinorWaveFunction>    & f1,
		  vector<SpinorBarWaveFunction> & a1,
		  vector<VectorWaveFunction>    & v1,
		  vector<VectorWaveFunction>    & v2,
		  bool me) const;
  
  /**
   * Matrix element for \f$f\bar{f}\to Z^0 \gamma\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  \f$Z^0\f$ wavefunction
   * @param v2  \f$\gamma\f$ wavefunction
   * @param me  Whether or not to calculate the matrix element for spin correlations
   */
  double ZGammaME(vector<SpinorWaveFunction>    & f1,
		  vector<SpinorBarWaveFunction> & a1,
		  vector<VectorWaveFunction>    & v1,
		  vector<VectorWaveFunction>    & v2,
		  bool me) const;

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
  static ClassDescription<MEPP2VGamma> initMEPP2VGamma;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEPP2VGamma & operator=(const MEPP2VGamma &) = delete;

private:

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
  AbstractFFVVertexPtr FFWvertex_;

  /**
   *   FFZVertex
   */
  AbstractFFVVertexPtr FFZvertex_;

  /**
   *  WWW Vertex
   */ 
  AbstractVVVVertexPtr WWWvertex_;
  //@}

  /**
   *  Processes
   */
  unsigned int process_;

  /**
   *  Allowed flavours of the incoming quarks
   */
  int maxflavour_;

  /**
   *  Treatment of the the boson masses
   */
  unsigned int massOption_;

  /**
   *  The matrix element
   */
  mutable ProductionMatrixElement me_;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MEPP2VGamma. */
template <>
struct BaseClassTrait<Herwig::MEPP2VGamma,1> {
  /** Typedef of the first base class of MEPP2VGamma. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEPP2VGamma class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEPP2VGamma>
  : public ClassTraitsBase<Herwig::MEPP2VGamma> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEPP2VGamma"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEPP2VGamma is implemented. It may also include several, space-separated,
   * libraries if the class MEPP2VGamma depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMEHadron.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEPP2VGamma_H */
