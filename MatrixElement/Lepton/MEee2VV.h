// -*- C++ -*-
#ifndef HERWIG_MEee2VV_H
#define HERWIG_MEee2VV_H
//
// This is the declaration of the MEee2VV class.
//

#include "Herwig/MatrixElement/HwMEBase.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "Herwig/MatrixElement/ProductionMatrixElement.h"

namespace Herwig {

using namespace ThePEG;

/**
 * The MEee2VV class implements the matrix elements for 
 * \f$e^+e^-\to W^+W^-/|^0Z^0\f$.
 *
 * @see \ref MEee2VVInterfaces "The interfaces"
 * defined for MEee2VV.
 */
class MEee2VV: public HwMEBase {

public:

  /**
   * The default constructor.
   */
  MEee2VV();

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
   * Used internally by generateKinematics, after calculating the
   * limits on cos(theta).
   */
  virtual double getCosTheta(double cthmin, double cthmax, const double r);

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
   * Matrix element for \f$f\bar{f}\to W^+W^-\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  Polarization vector for the 1st outgoing boson
   * @param v2  Polarization vector for the 2nd outgoing boson
   */
  double WWME(vector<SpinorWaveFunction>    & f1,
	      vector<SpinorBarWaveFunction> & a1,
	      vector<VectorWaveFunction>    & v1,
	      vector<VectorWaveFunction>    & v2) const;

  /**
   * Matrix element for \f$f\bar{f}\to Z^0Z^0\f$.
   * @param f1  Spinors for the incoming fermion
   * @param a1  Spinors for the incoming antifermion
   * @param v1  Polarization vector for the 1st outgoing boson
   * @param v2  Polarization vector for the 2nd outgoing boson
   */
  double ZZME(vector<SpinorWaveFunction>    & f1,
	      vector<SpinorBarWaveFunction> & a1,
	      vector<VectorWaveFunction>    & v1,
	      vector<VectorWaveFunction>    & v2) const;

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
  static ClassDescription<MEee2VV> initMEee2VV;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MEee2VV & operator=(const MEee2VV &) = delete;

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
   *  Treatment of the the W and Z masses
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
 *  base classes of MEee2VV. */
template <>
struct BaseClassTrait<Herwig::MEee2VV,1> {
  /** Typedef of the first base class of MEee2VV. */
  typedef Herwig::HwMEBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MEee2VV class and the shared object where it is defined. */
template <>
struct ClassTraits<Herwig::MEee2VV>
  : public ClassTraitsBase<Herwig::MEee2VV> {
  /** Return a platform-independent class name */
  static string className() { return "Herwig::MEee2VV"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MEee2VV is implemented. It may also include several, space-separated,
   * libraries if the class MEee2VV depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HwMELepton.so"; }
};

/** @endcond */

}

#endif /* HERWIG_MEee2VV_H */
